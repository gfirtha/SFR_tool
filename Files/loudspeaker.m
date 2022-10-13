classdef loudspeaker < handle
    %VIRTUAL_SOURCE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        source_index
        source_signal
        position
        orientation
        source_type
    end
    
    methods
        function obj = loudspeaker(idx, position, orientation, type)
            obj.source_index = idx;
            obj.position = position;
            obj.orientation = orientation;
            obj.source_signal = signal;
            obj.source_type = type;
        end
        
        function obj = set_input(obj,varargin)
            if length(varargin) == 1
                obj.source_signal.set_signal(varargin{1});
            else
                obj.source_signal.set_signal(varargin{1},varargin{2});
            end
        end
    end
end
