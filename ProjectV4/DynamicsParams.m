classdef DynamicsParams < handle
    properties
        lambda_front
        I_1
        I_2
        I_3
        I_4
        B_1
        B_2
        B_3
        B_4
        I_zz
        m_vehicle
        F_z1_static
        F_z2_static
        F_z3_static
        F_z4_static
    end

    methods
        function obj = DynamicsParams()
            obj.lambda_front = sym('lambda_front', 'real');
            obj.I_1          = sym('I_1', 'real');
            obj.I_2          = sym('I_2', 'real');
            obj.I_3          = sym('I_3', 'real');
            obj.I_4          = sym('I_4', 'real');

            obj.B_1 = sym('B_1', 'real');
            obj.B_2 = sym('B_2', 'real');
            obj.B_3 = sym('B_3', 'real');
            obj.B_4 = sym('B_4', 'real');

            obj.I_zz       = sym('I_zz', 'real');
            obj.m_vehicle  = sym('m_vehicle', 'real');

            obj.F_z1_static = sym('F_z1_static', 'real');
            obj.F_z2_static = sym('F_z2_static', 'real');
            obj.F_z3_static = sym('F_z3_static', 'real');
            obj.F_z4_static = sym('F_z4_static', 'real');
        end
    end
end
