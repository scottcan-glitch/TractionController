% % Inherits Car class. Add car presets here
% % 11/7: "sedan" preset added

classdef CarWithDefaults < Car
    methods
        function obj = CarWithDefaults(type)
            arguments
                type {mustBeMember(type,["sedan","custom"])} = "sedan"
            end

            obj@Car(); % call parent constructor

            if type == "sedan"
                obj = obj.setSedanDefaults();
            end
        end
    end

    methods (Access=private)
        function obj = setSedanDefaults(obj)

            % Pacejka
            obj.pacejka.long.B  = 10 * ones(4,1);
            obj.pacejka.long.C  = 1.9 * ones(4,1);
            obj.pacejka.long.E  = 0.97 * ones(4,1);
            obj.pacejka.long.mu = 0.9 * ones(4,1);

            obj.pacejka.lat.B   = 8 * ones(4,1);
            obj.pacejka.lat.C   = 1.3 * ones(4,1);
            obj.pacejka.lat.E   = 1.0 * ones(4,1);
            obj.pacejka.lat.mu  = 0.9 * ones(4,1);

            % Geometry
            obj.geometry.a       = 1.50;
            obj.geometry.b       = 1.20;
            obj.geometry.h       = 0.55;
            obj.geometry.t       = 1.60;
            obj.geometry.r_wheel = 0.32;

            % Dynamics
            obj.dynamics.lambda_front = 0.45;
            obj.dynamics.I_1 = 1.2;
            obj.dynamics.I_2 = 1.2;
            obj.dynamics.I_3 = 1.2;
            obj.dynamics.I_4 = 1.2;

            obj.dynamics.B_1 = 0.01;
            obj.dynamics.B_2 = 0.01;
            obj.dynamics.B_3 = 0.01;
            obj.dynamics.B_4 = 0.01;

            obj.dynamics.I_zz     = 2400;
            obj.dynamics.m_vehicle = 1000;  %change from 1500

            obj.dynamics.F_z1_static = 3500;
            obj.dynamics.F_z2_static = 3500;
            obj.dynamics.F_z3_static = 3500;
            obj.dynamics.F_z4_static = 3500;

            % Aero
            obj.aero.rho_air   = 1.225;
            obj.aero.C_drag    = 0.28;
            obj.aero.A_frontal = 2.2;
            obj.aero.A_side    = 2.5;

            % Limits (inputs, states)
            obj.limits.Input.Torque     = 200;
            obj.limits.Input.delta      = deg2rad(25);
            obj.limits.Input.delta_dot  = deg2rad(360);
            obj.limits.State.V_x        = 50;
            obj.limits.State.V_y        = 20;
            obj.limits.State.r          = deg2rad(180);
            obj.limits.State.omega_2    = (obj.limits.State.V_x/obj.geometry.r_wheel);
            obj.limits.State.omega_3    = (obj.limits.State.V_x/obj.geometry.r_wheel);

            obj.limits.Input.Lowerbound = -[obj.limits.Input.Torque; obj.limits.Input.Torque; obj.limits.Input.delta; obj.limits.Input.delta_dot];
            obj.limits.Input.Upperbound = [obj.limits.Input.Torque; obj.limits.Input.Torque; obj.limits.Input.delta; obj.limits.Input.delta_dot];
            obj.limits.State.Lowerbound = -[obj.limits.State.V_x; obj.limits.State.V_y; obj.limits.State.r; obj.limits.State.omega_2; obj.limits.State.omega_3];
            obj.limits.State.Upperbound = [obj.limits.State.V_x; obj.limits.State.V_y; obj.limits.State.r; obj.limits.State.omega_2; obj.limits.State.omega_3];

        end
    end
end
