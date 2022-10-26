classdef sound_scene_renderer < handle
    % sound_scene_renderer class contains one SFS_renderer instance for each
    % virtual source
    % TODO: make input, rendeder out and binauarl buses
    % and "wire" them up
    properties
        SFS_renderer
        directivity_tables
    end

    methods
        function obj = sound_scene_renderer(virtual_sources,loudspeaker,receiver, setup)

            % Get all required directivity characteristics
            N_fft = 2^nextpow2( min(setup.Block_size + size(setup.HRTF.Data.IR,3), 2*setup.Block_size) - 1 );
            cnt = 0;
            for n = 1 : length(virtual_sources)
                if isempty(get_dirtable_idx( obj.directivity_tables, virtual_sources{n}))
                    cnt = cnt + 1;
                    obj.directivity_tables{cnt} = directivity_table(virtual_sources{n}.source_type, N_fft, setup.Input_stream.SampleRate);
                end
            end

            for n = 1 : length(virtual_sources)
                idx = get_dirtable_idx(obj.directivity_tables,virtual_sources{n});
                switch virtual_sources{n}.renderer_type
                    case 'VBAP'
                        obj.SFS_renderer{n} = vbap_renderer(virtual_sources{n}, loudspeaker);
                    case 'DBAP'
                        obj.SFS_renderer{n} = dbap_renderer(virtual_sources{n}, loudspeaker);
                    case 'WFS'
                        obj.SFS_renderer{n} = wfs_renderer(virtual_sources{n}, loudspeaker, setup.Input_stream.SampleRate,obj.directivity_tables{idx}, setup.Renderer_setup.Antialiasing);
                    case 'TD_stereo'
                        obj.SFS_renderer{n} = time_delay_renderer(virtual_sources{n}, loudspeaker, setup.Input_stream.SampleRate);
                    case 'CTC'
                        obj.SFS_renderer{n} = ctc_renderer(virtual_sources{n}, loudspeaker, receiver, setup.Input_stream.SampleRate, setup.Renderer_setup.Plant_model, setup.Renderer_setup.VS_model, setup.Renderer_setup.HRTF_database,setup.Renderer_setup.N_filt);
                end
            end
        end

        function update_SFS_renderers(obj, type)
            for n = 1 : length(obj.SFS_renderer)
                obj.SFS_renderer{n}.update_renderer(type);
            end
        end

        function render(obj, input)
            for m = 1 : length(obj.SFS_renderer)
                obj.SFS_renderer{m}.virtual_source.source_signal.set_signal(input(:,m));
                obj.SFS_renderer{m}.render;
                for n = 1 : length(obj.SFS_renderer{m}.secondary_source_distribution)
                    obj.SFS_renderer{m}.secondary_source_distribution{n}.output_signal.add_signals(...
                        obj.SFS_renderer{m}.output_signal{n});
                end
            end
        end
    end
end