function [Slip, FrontWheels, wc] = makeKinematics(x,x_dot,u,u_dot,car)

% % Returns symbolic Slip structure
% % Steer angle, wheel angular velocity, and vehicle frame velocities are
% % known, slips can be resolved symbolically via kinematics
% return this also?: WheelCenterVelocity, WheelCenterAcceleration,

% % EXTRACT

    % Extract states for slip calcs
    omega_2 = x.omega_2;
    omega_3 = x.omega_3;
    v_x = x.v_x;
    v_y = x.v_y;
    r = x.r;
    
    % Extract StateDerivs
    v_dot_x = x_dot.v_dot_x;
    v_dot_y = x_dot.v_dot_y;
    r_dot = x_dot.r_dot;
    
    % Extract car parameters to define wheel locations (&radius)
    wc = struct(); % return structure
    wc.x = [car.geometry.b,-car.geometry.a,-car.geometry.a,car.geometry.b];
    wc.y = [car.geometry.t/2,car.geometry.t/2,-car.geometry.t/2,-car.geometry.t/2];

    wc_1_x = wc.x(1);
    wc_1_y = wc.y(1);
    wc_2_x = wc.x(2);
    wc_2_y = wc.y(2);
    wc_3_x = wc.x(3);
    wc_3_y = wc.y(3);
    wc_4_x = wc.x(4);
    wc_4_y = wc.y(4);
    r_wheel = car.geometry.r_wheel;

    % Extract Input
    delta = u.delta;

% % END

% Need to take time deriv of WheelFrameWheelVelocity, use these helper
% functions for product/chain rule, includes steer rate terms d(delta)/dt
% Ai, Bi are the terms used in wheel center velocity. ex:
% wc_1_vx = A1*cos(delta) + B1*sin(delta)

    A1 = v_x - wc_1_y*r;   B1 = v_y + wc_1_x*r;
    A2 = v_x - wc_2_y*r;   B2 = v_y + wc_2_x*r;
    A3 = v_x - wc_3_y*r;   B3 = v_y + wc_3_x*r;
    A4 = v_x - wc_4_y*r;   B4 = v_y + wc_4_x*r;
    A1_dot = v_dot_x - wc_1_y*r_dot - wc_1_x*r^2;
    B1_dot = v_dot_y + wc_1_x*r_dot - wc_1_y*r^2;
    % A2_dot = v_dot_x - wc_2_y*r_dot - wc_2_x*r^2;
    % B2_dot = v_dot_y + wc_2_x*r_dot - wc_2_y*r^2;
    % A3_dot = v_dot_x - wc_3_y*r_dot - wc_3_x*r^2;
    % B3_dot = v_dot_y + wc_3_x*r_dot - wc_3_y*r^2;
    A4_dot = v_dot_x - wc_4_y*r_dot - wc_4_x*r^2;
    B4_dot = v_dot_y + wc_4_x*r_dot - wc_4_y*r^2;
% WheelFrameWheelVelocity
wc_1_vx = A1*cos(delta) + B1*sin(delta);
wc_1_vy = B1*cos(delta) - A1*sin(delta);
wc_2_vx = A2;
wc_2_vy = B2;
wc_3_vx = A3;
wc_3_vy = B3;
wc_4_vx = A4*cos(delta) + B4*sin(delta);
wc_4_vy = B4*cos(delta) - A4*sin(delta);

% WheelFrameWheelAcceleration
    % Only front wheel x direction accel is needed.
    % Others are here for fun & commented out

    % Front (steered) wheels
wc_1_ax = A1_dot*cos(delta) - A1*sin(delta)*u_dot.delta_dot + B1_dot*sin(delta) + B1*cos(delta)*u_dot.delta_dot;
    % wc_1_ay = B1_dot*cos(delta) - B1*x_dot.delta_dot*sin(delta) - (A1_dot*sin(delta) + A1*x_dot.delta_dot*cos(delta));
wc_4_ax = A4_dot*cos(delta) - A4*u_dot.delta_dot*sin(delta) + B4_dot*sin(delta) + B4*u_dot.delta_dot*cos(delta);
    % wc_4_ay = B4_dot*cos(delta) - B4*x_dot.delta_dot*sin(delta) - (A4_dot*sin(delta) + A4*x_dot.delta_dot*cos(delta));

    % Rear (unsteered) wheels - 
    % wc_2_ax = A2_dot;
    % wc_2_ay = B2_dot;
    % wc_3_ax = A3_dot;
    % wc_3_ay = B3_dot;


% ZeroLongSlipFront
omega_1 = wc_1_vx / r_wheel;
omega_4 = wc_4_vx / r_wheel;
omega_dot_1 = wc_1_ax / r_wheel;
omega_dot_4 = wc_4_ax / r_wheel;

FrontWheels = struct();
FrontWheels.omega_1 = omega_1;
FrontWheels.omega_4 = omega_4;
FrontWheels.omega_dot_1 = omega_dot_1;
FrontWheels.omega_dot_4 = omega_dot_4;

% Be careful of divide by zero here, might need max("",1e-6) in den
% LongitudinalSlip

epsVal = 1e-6; %avoids divide by zero, & doesnt break symbolic solvers
% wheel center velocity in x must be positive

% % REMOVED ABS(WC_I_VX)?????????????????
lambda_1 = (omega_1*r_wheel - wc_1_vx) ./ (wc_1_vx + epsVal);
lambda_2 = (omega_2*r_wheel - wc_2_vx) ./ (wc_2_vx + epsVal);
lambda_3 = (omega_3*r_wheel - wc_3_vx) ./ (wc_3_vx + epsVal);
lambda_4 = (omega_4*r_wheel - wc_4_vx) ./ (wc_4_vx + epsVal);

% LateralSlipAngle
alpha_1 = atan(wc_1_vy/wc_1_vx);
alpha_2 = atan(wc_2_vy/wc_2_vx);
alpha_3 = atan(wc_3_vy/wc_3_vx);
alpha_4 = atan(wc_4_vy/wc_4_vx);


Slip = struct();

Slip.long = [lambda_1, lambda_2, lambda_3, lambda_4];
Slip.lat = [alpha_1, alpha_2, alpha_3, alpha_4];


end