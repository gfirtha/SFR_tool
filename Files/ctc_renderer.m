classdef ctc_renderer < handle
    %CTC_RENDERER Summary of this class goes here
    %   Detailed explanation goes here
    %
    %
    %   TODO: implement other models for plant_model and
    %   virtual_source_model
    %
    %
    properties
        virtual_source
        secondary_source_distribution
        receiver
        output_signal
        fs
        plant_model
        virtual_source_model    % Virtual source to receiver propagation in frequency domain
        hrtf_database
        hrtf_2d_database
        inv_plant_mx_f    % Inverse plant matrix in frequency domain
        virtual_source_coefficients
        N_filt
    end
    properties (SetAccess = protected)
        CTC_filters
    end

    methods
        function obj = ctc_renderer(varargin)
            obj.virtual_source = varargin{1};
            obj.secondary_source_distribution = varargin{2};
            obj.receiver = varargin{3};
            obj.fs = varargin{4};
            for n = 1 : length(obj.secondary_source_distribution)
                obj.output_signal{n} = signal;
            end
            obj.plant_model = varargin{5};
            obj.virtual_source_model = varargin{6};
            if (strcmp(obj.plant_model,'HRTF') || strcmp(obj.virtual_source_model,'HRTF'))
                obj.hrtf_database = varargin{7};
                obj.N_filt = size(obj.hrtf_database.Data.IR,3);

                R_measurement = mean(obj.hrtf_database.SourcePosition(:,3));
                theta_measurement = obj.hrtf_database.SourcePosition(:,1:2);
                ixs = find(theta_measurement(:,2) == 0);
                theta_measurement = theta_measurement(ixs,1)*pi/180;
                hrtf_measured = obj.hrtf_database.Data.IR(ixs,:,:);
                [theta_measurement,ix] = sort(theta_measurement);
                hrtf_measured = hrtf_measured(ix,:,:);
              %  theta_measurement = theta_measurement(1:16:end);
              %  hrtf_measured = hrtf_measured(1:16:end,:,:);
                obj.hrtf_2d_database = struct('R',R_measurement,'theta',theta_measurement, 'spectrum',fft(hrtf_measured,[],3));

            else
                obj.N_filt = varargin{8};
            end

            obj.update_plant_mx;
            obj.update_vs_model;
            driving_signal = obj.get_driving_filter;
            for n = 1 : length(obj.secondary_source_distribution)
                obj.CTC_filters{n}  = OLS_convolver(driving_signal (:,n), length(obj.virtual_source.source_signal.time_series));
            end

        end


        function obj = update_plant_mx(obj)
            switch obj.plant_model
                case 'HRTF'
                    xs = cell2mat(cellfun( @(x) x.position,    obj.secondary_source_distribution, 'UniformOutput', false)');
                    plant_mx_t = get_hrtfs( xs, obj.receiver.position, obj.receiver.orientation, obj.hrtf_database, obj.hrtf_2d_database );
                    plant_mx_f = fft(plant_mx_t,[],3);
                    f = reshape((0:obj.N_filt-1)'/obj.N_filt*obj.fs, [1,1,obj.N_filt] ) ;
                case 'point_source'
                    xs = cell2mat(cellfun( @(x) x.position,    obj.secondary_source_distribution, 'UniformOutput', false)');
                    r_head = 0.1;
                    x_ear = bsxfun( @plus, obj.receiver.position', fliplr((obj.receiver.orientation'*[1,-1]*r_head)'));
                    
                    vs_pos=obj.virtual_source.position;
                    receiver_pos=obj.receiver.position;
                    receiver_dir=rad2deg(obj.receiver.orientation);
                    r_head = 0.1;
                    L_ear_pos = receiver_pos + abs(deg2rad(receiver_dir+90))*r_head;
                    R_ear_pos = receiver_pos + abs(deg2rad(receiver_dir-90))*r_head;
                    
                    Rmx = zeros(size(x_ear,1), size(xs,1));
                    for n = 1 : size(xs,1)
                        v_sr = bsxfun(@minus, xs(n,:), x_ear);
                        Rmx(:,n) = sqrt(sum(v_sr.^2,2));
                    end
                    f = reshape((0:obj.N_filt-1)'/obj.N_filt*obj.fs, [1,1,obj.N_filt] ) ;
                    plant_mx_f = 1/(4*pi)*bsxfun( @times, exp( -1i*2*pi*bsxfun( @times, f, Rmx/340  ) ), 1./Rmx);

% =======
%                     vs_pos=obj.virtual_source.position;
%                     receiver_pos=obj.receiver.position;
%                     receiver_dir=rad2deg(obj.receiver.orientation);
%                     r_head = 0.1;
%                     L_ear_pos = receiver_pos + abs(deg2rad(receiver_dir+90))*r_head;
%                     R_ear_pos = receiver_pos + abs(deg2rad(receiver_dir-90))*r_head;
%                     dist_table = [ norm(L_ear_pos - vs_pos) ; norm(R_ear_pos - vs_pos)];
%                     dist_center = norm(receiver_pos-vs_pos)
%                     c=environment.c;
%                     M = 4096;
%                     f = (0:M/2-1)/M*fs;
%                     for fi = 1 : length(f)
%                         plant_mx_f(:,:,fi) = exp(-1i*2*pi*f(fi)*dist_table/c)./dist_table.*dist_center./exp(-1i*2*pi*f(fi)*dist_center/c);
%                         obj.inv_plant_mx_f(:,:,fi) = pinv(plant_mx_f(:,:,fi));
%                     end
% >>>>>>> origin/kissdani
            end

            obj.inv_plant_mx_f = zeros(size(plant_mx_f));
            for n = 1 : size(plant_mx_f,3)
                X = squeeze(plant_mx_f(:,:,n));
                lambda = 1e-5;
                obj.inv_plant_mx_f(:,:,n) = inv(X.'*X + lambda*eye(size(X)))*X.';
            end
        end

        function obj = update_vs_model(obj)
            switch obj.virtual_source_model
                case 'HRTF'
                    obj.virtual_source_coefficients = fft(get_hrtfs( obj.virtual_source.position, obj.receiver.position, obj.receiver.orientation, obj.hrtf_database, obj.hrtf_2d_database  ),[],2);
                case 'point_source'
                    xs = obj.virtual_source.position;
                    r_head = 0.1;
                    x_ear = bsxfun( @plus, obj.receiver.position', fliplr((obj.receiver.orientation'*[1,-1]*r_head)'));
                    R = sqrt(sum( (bsxfun( @plus, x_ear, -xs)).^2,2));
                    f = (0:obj.N_filt-1)'/obj.N_filt*obj.fs ;
                    obj.virtual_source_coefficients = 1/(4*pi)* bsxfun(@times, exp( -1i*2*pi* bsxfun(@times, f', R/340 )), 1./R );
            end
        end

        function obj = update_renderer(obj,type)
            switch type
                case 'receiver_moved'
                    obj.update_plant_mx;
                    obj.update_vs_model;
                case 'receiver_rotated'
                    obj.update_plant_mx;
                    obj.update_vs_model;
                case 'loudspeaker_moved'
                    obj.update_plant_mx;
                case 'loudspeaker_rotated'
                    obj.update_plant_mx;
                case 'virtual_source_moved'
                    obj.update_vs_model;
                case 'virtual_source_rotated'
                    obj.update_vs_model;
            end
            driving_signal = obj.get_driving_filter;
            for n = 1 : length(obj.secondary_source_distribution)
                obj.CTC_filters{n}.update_coefficients(driving_signal (:,n));
            end

        end

        function driving_filter = get_driving_filter(obj)
            N = size(obj.inv_plant_mx_f,3);
            driving_function = zeros(N, size(obj.inv_plant_mx_f,2));
            for n = 1 : N
                driving_function(n,:) = squeeze(obj.inv_plant_mx_f(:,:,n))*obj.virtual_source_coefficients(:,n);
            end
            driving_filter = fftshift(ifft(driving_function,N,1,'symmetric'),1).*tukeywin(N,0.1);
        end

        function render(obj)
            for n = 1 : length(obj.secondary_source_distribution)
                obj.output_signal{n}.set_signal( obj.CTC_filters{n}.convolve( obj.virtual_source.source_signal.time_series ) );
            end
        end

    end
end

