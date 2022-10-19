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
                    xs = cell2mat(cellfun( @(x) x.position,    obj.secondary_source_distribution, 'UniformOutput', false)');
                    plant_mx_t = get_hrtfs( xs, obj.receiver.position, obj.receiver.orientation, obj.hrtf_database );
                    plant_mx_f = fft(plant_mx_t,[],3);
                    obj.inv_plant_mx_f = zeros(size(plant_mx_f));
                    for n = 1 : size(plant_mx_t,3)
                        obj.inv_plant_mx_f(:,:,n) = pinv(squeeze(plant_mx_f(:,:,n)));
                    end
                case 'point_source'
            end

        end

        function obj = update_vs_model(obj)
            switch obj.virtual_source_model
                case 'HRTF'
                    obj.virtual_source_coefficients = fft(get_hrtfs( obj.virtual_source.position, obj.receiver.position, obj.receiver.orientation, obj.hrtf_database ),[],2);
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

