classdef Limits < handle
    properties
        % Input
        Input
        Torque   %for left and right rear wheel motors
        delta
        delta_dot

        % State
        State
        V_x
        V_y
        r
        omega_2
        omega_3

        % vectors
        Lowerbound
        Upperbound

    end

    methods
        function obj = Limits()
            obj.Input.Torque  = sym('Torque', 'real');
            obj.Input.delta  = sym('delta', 'real');
            
            obj.State.V_x = sym('V_x', 'real');
            obj.State.V_y = sym('V_y', 'real');
            obj.State.r = sym('r', 'real');
            obj.State.omega_2= sym('omega_2', 'real');
            obj.State.omega_3 = sym('omega_3', 'real');

            obj.Input.Lowerbound = -[obj.Input.Torque; obj.Input.Torque; obj.Input.delta];
            obj.Input.Upperbound = [obj.Input.Torque; obj.Input.Torque; obj.Input.delta];
            obj.State.Lowerbound = -[obj.State.V_x; obj.State.V_y; obj.State.r; obj.State.omega_2; obj.State.omega_3];
            obj.State.Upperbound = [obj.State.V_x; obj.State.V_y; obj.State.r; obj.State.omega_2; obj.State.omega_3];

        end
    end
end
