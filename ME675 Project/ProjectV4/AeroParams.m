classdef AeroParams < handle
    properties
        rho_air
        C_drag
        A_frontal
        A_side
    end

    methods
        function obj = AeroParams()
            obj.rho_air  = sym('rho_air', 'real');
            obj.C_drag   = sym('C_drag', 'real');
            obj.A_frontal = sym('A_frontal', 'real');
            obj.A_side   = sym('A_side', 'real');
        end
    end
end
