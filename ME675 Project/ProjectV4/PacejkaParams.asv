classdef PacejkaParams < handle
    properties
        % Longitudinal
        long
        B_long
        C_long
        E_long
        mu_long

        % Lateral
        lat
        B_lat
        C_lat
        E_lat
        mu_lat

    end

    methods
        function obj = PacejkaParams()
            obj.long.B  = sym('B_long',  [4 1], 'real');
            obj.long.C  = sym('C_long',  [4 1], 'real');
            obj.long.E  = sym('E_long',  [4 1], 'real');
            obj.long.mu = sym('mu_long', [4 1], 'real');

            obj.lat.B   = sym('B_lat',   [4 1], 'real');
            obj.lat.C   = sym('C_lat',   [4 1], 'real');
            obj.lat.E   = sym('E_lat',   [4 1], 'real');
            obj.lat.mu  = sym('mu_lat',  [4 1], 'real');
        end
    end
end
