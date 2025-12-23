classdef GeometryParams < handle
    properties
        a
        b
        h
        t
        r_wheel
    end

    methods
        function obj = GeometryParams()
            obj.a       = sym('a', 'real');
            obj.b       = sym('b', 'real');
            obj.h       = sym('h', 'real');
            obj.t       = sym('t', 'real');
            obj.r_wheel = sym('r_wheel', 'real');
        end
    end
end
