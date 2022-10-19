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
        inv_plant_mx_f    % Inverse plant matrix in frequency domain
        virtual_source_coefficients
    end
    properties (SetAccess = protected)
        CTC_filters
    end

    methods
        function obj = ctc_renderer(virtual_source,SSD,receiver, fs,mode1,mode2,hrtf_in)
            obj.virtual_source = virtual_source;
            obj.secondary_source_distribution = SSD;
            obj.receiver = receiver;
            obj.fs = fs;
            for n = 1 : length(obj.secondary_source_distribution)
                obj.output_signal{n} = signal;
            end
            obj.plant_model = mode1;
            obj.virtual_source_model = mode2;
            obj.hrtf_database = hrtf_in;
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
                    x0 = cell2mat(cellfun( @(x) x.position,    obj.secondary_source_distribution, 'UniformOutput', false)');
                    v_vec = bsxfun( @minus, x0, obj.receiver.position);
                    v_vec = bsxfun(@times, v_vec, 1./sqrt(sum(v_vec.^2,2)));
                    [theta_1,~] = cart2pol(v_vec(:,1), v_vec(:,2));
                    [theta_2,~] = cart2pol(obj.receiver.orientation(1),obj.receiver.orientation(2));
                    for n = 1 : length(theta_1)
                        [~,inds(n)] = min( sum(  (bsxfun(@minus, obj.hrtf_database.SourcePosition(:,[1,2]) ...
                            ,mod([theta_1(n) - theta_2,0]*180/pi, 360))).^2,2  ) );
                    end
                    % row index: ear,  column index: loudspeaker
                    plant_mx_t = zeros(2, length(obj.secondary_source_distribution),size(obj.hrtf_database.Data.IR,3));
                    plant_mx_t(1,:,:) = squeeze(obj.hrtf_database.Data.IR(inds,1,:));
                    plant_mx_t(2,:,:) = squeeze(obj.hrtf_database.Data.IR(inds,2,:));
                    plant_mx_f = fft(plant_mx_t,[],3);
                    obj.inv_plant_mx_f = zeros(size(plant_mx_f));
                    for n = 1 : size(plant_mx_t,3)
                        obj.inv_plant_mx_f(:,:,n) = pinv(squeeze(plant_mx_f(:,:,n)));
                    end
                case 'point_source'
                    obj;
            end

        end

        function obj = update_vs_model(obj)

            switch obj.virtual_source_model
                case 'HRTF'
                    xs_vec = obj.virtual_source.position-obj.receiver.position;
                    xs_vec = bsxfun(@times, xs_vec, 1./sqrt(sum(xs_vec.^2,2)));
                    [theta_1,~] = cart2pol(xs_vec(:,1), xs_vec(:,2));
                    [theta_2,~] = cart2pol(obj.receiver.orientation(1),obj.receiver.orientation(2));
                    [~,ind] = min( sum(  (bsxfun(@minus, obj.hrtf_database.SourcePosition(:,[1,2]) ...
                        ,mod([theta_1 - theta_2,0]*180/pi, 360))).^2,2  ) );
                    obj.virtual_source_coefficients = fft([squeeze(obj.hrtf_database.Data.IR(ind,1,:)),...
                        squeeze(obj.hrtf_database.Data.IR(ind,2,:))],[],1).';
                case 'point_source'
            end
        end

        function obj = update_renderer(obj)
            obj.update_plant_mx;
            obj.update_vs_model;
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
            driving_filter = fftshift(ifft(driving_function,N,1,'symmetric'),1).*tukeywin(N,0.35);
        end

        function render(obj)
            for n = 1 : length(obj.secondary_source_distribution)
                obj.output_signal{n}.set_signal( obj.CTC_filters{n}.convolve( obj.virtual_source.source_signal.time_series ) );
            end
        end

    end
end

