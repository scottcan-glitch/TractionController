function [omega_dot_2, omega_dot_3] = makeDriveTrain(x, u, car, WheelFrameForces)

% % EXTRACT
    % Extract from WheelFrameForces
    F_2_long = WheelFrameForces.long(2);
    F_3_long = WheelFrameForces.long(3);
    
    % Extract Input
    T_m2 = u.T_m2;
    T_m3 = u.T_m3;
    
    % Extract State
    omega_2 = x.omega_2;
    omega_3 = x.omega_3;
    
    % Extract from parameters
    B_2 = car.dynamics.B_2;
    B_3 = car.dynamics.B_3;
    I_2 = car.dynamics.I_2;
    I_3 = car.dynamics.I_3;
    r_wheel = car.geometry.r_wheel;

% % END

% PowerTrainEquations
omega_dot_2 = (T_m2 - F_2_long*r_wheel - B_2*omega_2)/I_2;
omega_dot_3 = (T_m3 - F_3_long*r_wheel - B_3*omega_3)/I_3;


F_1_long = WheelFrameForces.long(1);
F_4_long = WheelFrameForces.long(4);

omega_dot_1 = (F_1_long*r_wheel - B_1*omega_2)/I_1;
omega_dot_4 = (F_4_long*r_wheel - B_4*omega_4)/I_4;
end