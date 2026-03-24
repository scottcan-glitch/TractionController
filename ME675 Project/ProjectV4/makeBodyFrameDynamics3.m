function [v_dot_x, v_dot_y, r_dot] = makeBodyFrameDynamics3(car, State, StateDeriv, Input, WheelFrameForces, wc)

% % Function uses this model:
% (v_dot_x - v_y*r)*m_vehicle == F_x;
% (v_dot_y + v_x*r)*m_vehicle == F_y;
% r_dot*I_zz == M_z;

    % % EXTRACT
    % Extract States
    v_x = State.v_x;
    v_y = State.v_y;
    r = State.r;
    v_dot_x = StateDeriv.v_dot_x;
    v_dot_y = StateDeriv.v_dot_y;
    r_dot = StateDeriv.r_dot;
    
    % % Extract Parameters
    m_vehicle = car.dynamics.m_vehicle;
    I_zz = car.dynamics.I_zz;
    delta = Input.delta;
    
    % wheel loc for moments (from wc struct, kinematics)
    wc_1_x = wc.x(1);
    wc_1_y = wc.y(1);
    wc_2_x = wc.x(2);
    wc_2_y = wc.y(2);
    wc_3_x = wc.x(3);
    wc_3_y = wc.y(3);
    wc_4_x = wc.x(4);
    wc_4_y = wc.y(4);

    % Extract aero car parameters
    C_drag = car.aero.C_drag;
    rho_air = car.aero.rho_air;
    A_frontal = car.aero.A_frontal;
    A_side = car.aero.A_side;

% % END


% Helper expressions
F_x1 = WheelFrameForces.long(1)*cos(delta) - WheelFrameForces.lat(1)*sin(delta);
F_y1 = WheelFrameForces.lat(1)*cos(delta) + WheelFrameForces.long(1)*sin(delta);
F_x2 = WheelFrameForces.long(2);
F_y2 = WheelFrameForces.lat(2);
F_x3 = WheelFrameForces.long(3);
F_y3 = WheelFrameForces.lat(3);
F_x4 = WheelFrameForces.long(4)*cos(delta) - WheelFrameForces.lat(4)*sin(delta);
F_y4 = WheelFrameForces.lat(4)*cos(delta) + WheelFrameForces.long(4)*sin(delta);
F_drag_x = 0.5*C_drag*rho_air*v_x^2*A_frontal;
F_drag_y = 0.5*C_drag*rho_air*v_y^2*A_side;

disp('sums:')
% Sums
F_x = F_x1 + F_x2 + F_x3 + F_x4 - F_drag_x
F_y = F_y1 + F_y2 + F_y3 + F_y4 - F_drag_y
M_z = wc_1_x * F_y1 - wc_1_y * F_x1 + wc_2_x * F_y2 - ...
    wc_2_y * F_x2 + wc_3_x * F_y3 - wc_3_y * F_x3 + ...
    wc_4_x * F_y4 - wc_4_y * F_x4

symvar(F_x)
symvar(F_y)
symvar(M_z)

% Newtons 3DOF, stored for linear solve w/ z
v_dot_x = (F_x/m_vehicle) + v_y*r
v_dot_y = (F_y/m_vehicle) - v_x*r
r_dot = M_z/I_zz

end