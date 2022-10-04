classdef environment
    %ENVIRONMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        rho
        c
    end
    
    methods
        function obj = environment
            %ENVIRONMENT Construct an instance of this class
            %   if the simulaiton would take place under different
            %   cirumstances, change this: eg.:
            %   underwater: rho = 1000 ; c = 1450
            obj.rho = 1.2;
            obj.c = 343.1;
        end
        
    end
end

