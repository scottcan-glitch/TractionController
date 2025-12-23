classdef Kart < handle
    properties
        pacejka   % PacejkaParams
        geometry  % GeometryParams
        dynamics  % DynamicsParams
        aero      % AeroParams
        limits    % Limits (state and input constraints)
    end

    methods
        function obj = Kart()
            obj.pacejka   = PacejkaParams();
            obj.geometry  = GeometryParams();
            obj.dynamics  = DynamicsParams();
            obj.aero      = AeroParams();
            obj.limits    = Limits_3input();
        end
    end
end
