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
        r_head = 0.1 % radius of head for rigid sphere and point source model
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
                    freq = reshape((0:obj.N_filt-1)'/obj.N_filt*obj.fs, [1,1,obj.N_filt] ) ;
                case 'point_source'
                    xs = cell2mat(cellfun( @(x) x.position,    obj.secondary_source_distribution, 'UniformOutput', false)');
                    %r_head = 0.1;
                    x_ear = bsxfun( @plus, obj.receiver.position', fliplr((obj.receiver.orientation'*[1,-1]*r_head)'));

                    Rmx = zeros(size(x_ear,1), size(xs,1));
                    for n = 1 : size(xs,1)
                        v_sr = bsxfun(@minus, xs(n,:), x_ear);
                        Rmx(:,n) = sqrt(sum(v_sr.^2,2));
                    end

                    freq = reshape((0:obj.N_filt-1)'/obj.N_filt*obj.fs, [1,1,obj.N_filt] ) ;
                    plant_mx_f = 1/(4*pi)*bsxfun( @times, exp( -1i*2*pi*bsxfun( @times, freq, Rmx/340  ) ), 1./Rmx);
                case 'rigid_sphere'
                    % A comparison of the performance of HRTF models in inverse filter design for Crosstalk Cancellation
                    % 2.2 (5)
%                     theta = obj.receiver.postion
%                     cLear = [];
%                     cRear = [];
%                     posSpeak{1}= obj.secondary_source_distribution{1, 1}.position;
%                     posSpeak{2}= obj.secondary_source_distribution{1, 2}.position;
                    
                    xs = cell2mat(cellfun( @(x) x.position,    obj.secondary_source_distribution, 'UniformOutput', false)');
                    %r_head = 0.1;
                    x_ear = bsxfun( @plus, obj.receiver.position', fliplr((obj.receiver.orientation'*[1,-1]*obj.r_head)'));

                    Rmx = zeros(size(x_ear,1), size(xs,1));                     
                    theta_mx = zeros(size(x_ear,1), size(xs,1));
                    for n = 1 : size(xs,1)
                        v_sr = bsxfun(@minus, xs(n,:), x_ear);
                        Rmx(:,n) = sqrt(sum(v_sr.^2,2));

                        v_ls = bsxfun(@minus, xs(n,:), obj.receiver.position);
                        v_ls/norm(v_ls);
                        v_ears = bsxfun(@minus, x_ear, obj.receiver.position);
                        v_ears = bsxfun(@times, v_ears, 1./sqrt(sum(v_ears.^2,2)));

                        theta_mx(:,n) = acos(v_ls*v_ears.');
                    end
                    c = 343.1;
                    freq = reshape((0:obj.N_filt-1)'/obj.N_filt*obj.fs, [1,1,obj.N_filt] ) ;
                    k = 2*pi*freq / c;
                    Norder = 50;
                     
                    plant_mx = zeros(2,size(xs,1),length(freq));
                    sign_mx = [-ones(size(xs,1));size(xs,1)];
                    for fi = 2 : length(freq)
                        A0 = Rmx./(k(fi)*obj.r_head^2.*exp(-1i*k(fi).*Rmx));
                        for n = 1 : Norder
                            plant_mx = plant_mx + A0*(2*n+1)*getSphH( n, 2, k(fi)*Rmx ).*sign_mx^n.*legendreP(n,2*cos(theta_mx))./getDifSphH( n, 2, k(fi)*Rmx );
                        end
                   end
            end

%             obj.inv_plant_mx_f = zeros(size(plant_mx_f));
%             for n = 1 : size(plant_mx_f,3)
%                 X = squeeze(plant_mx_f(:,:,n));
%                 lambda = 1e-5;
%                 obj.inv_plant_mx_f(:,:,n) = pinv(X.'*X + lambda*eye(size(X)))*X.';
%             end
%             obj.inv_plant_mx_f(:,:,squeeze(f)>20e3) = 0;
%             obj.inv_plant_mx_f = 10*obj.inv_plant_mx_f / max(max(max(obj.inv_plant_mx_f)));
        end

        function obj = update_vs_model(obj)
            switch obj.virtual_source_model
                case 'HRTF'
                    obj.virtual_source_coefficients = fft(get_hrtfs( obj.virtual_source.position, obj.receiver.position, obj.receiver.orientation, obj.hrtf_database, obj.hrtf_2d_database  ),[],2);
                case 'point_source'
                    xs = obj.virtual_source.position;
                    %r_head = 0.1;
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

