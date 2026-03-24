function [Slip, FrontWheels] = makeKinematics(x,x_dot,u,u_dot,P)

% % Returns symbolic Slip structure
% % Steer angle, wheel angular velocity, and vehicle frame velocities are
% % known, slips can be resolved symbolically via kinematics
% return this also?: WheelCenterVelocity, WheelCenterAcceleration,

% Extract states for slip calcs
omega_2 = x.omega_2;
omega_3 = x.omega_3;

% BodyFrameWheelLocations
wc_1_x = P.b;
wc_1_y = P.t/2;
wc_2_x = -P.a;
wc_2_y = P.t/2;
wc_3_x = -P.a;
wc_3_y = -P.t/2;
wc_4_x = P.b;
wc_4_y = -P.t/2;

% Need to take time deriv of WheelFrameWheelVelocity, use these helper
% functions for product/chain rule, includes steer rate terms d(delta)/dt
% Ai, Bi are the terms used in wheel center velocity. ex:
% wc_1_vx = A1*cos(delta) + B1*sin(delta)

    A1 = x.v_x - wc_1_y*x.r;   B1 = x.v_y + wc_1_x*x.r;
    A2 = x.v_x - wc_2_y*x.r;   B2 = x.v_y + wc_2_x*x.r;
    A3 = x.v_x - wc_3_y*x.r;   B3 = x.v_y + wc_3_x*x.r;
    A4 = x.v_x - wc_4_y*x.r;   B4 = x.v_y + wc_4_x*x.r;
    A1_dot = x_dot.v_dot_x - wc_1_y*x_dot.r_dot - wc_1_x*x.r^2;
    B1_dot = x_dot.v_dot_y + wc_1_x*x_dot.r_dot - wc_1_y*x.r^2;
    A2_dot = x_dot.v_dot_x - wc_2_y*x_dot.r_dot - wc_2_x*x.r^2;
    B2_dot = x_dot.v_dot_y + wc_2_x*x_dot.r_dot - wc_2_y*x.r^2;
    A3_dot = x_dot.v_dot_x - wc_3_y*x_dot.r_dot - wc_3_x*x.r^2;
    B3_dot = x_dot.v_dot_y + wc_3_x*x_dot.r_dot - wc_3_y*x.r^2;
    A4_dot = x_dot.v_dot_x - wc_4_y*x_dot.r_dot - wc_4_x*x.r^2;
    B4_dot = x_dot.v_dot_y + wc_4_x*x_dot.r_dot - wc_4_y*x.r^2;

% WheelFrameWheelVelocity
wc_1_vx = A1*cos(u.delta) + B1*sin(u.delta);
wc_1_vy = B1*cos(u.delta) - A1*sin(u.delta);
wc_2_vx = A2;
wc_2_vy = B2;
wc_3_vx = A3;
wc_3_vy = B3;
wc_4_vx = A4*cos(u.delta) + B4*sin(u.delta);
wc_4_vy = B4*cos(u.delta) - A4*sin(u.delta);

% WheelFrameWheelAcceleration
    % Only front wheel x direction accel is needed.
    % Others are here for fun & commented out

    % Front (steered) wheels
wc_1_ax = A1_dot*cos(u.delta) - A1*sin(u.delta)*u_dot.delta_dot + B1_dot*sin(u.delta) + B1*cos(u.delta)*u_dot.delta_dot;
    % wc_1_ay = B1_dot*cos(u.delta) - B1*x_dot.delta_dot*sin(u.delta) - (A1_dot*sin(u.delta) + A1*x_dot.delta_dot*cos(u.delta));
wc_4_ax = A4_dot*cos(u.delta) - A4*u_dot.delta_dot*sin(u.delta) + B4_dot*sin(u.delta) + B4*u_dot.delta_dot*cos(u.delta);
    % wc_4_ay = B4_dot*cos(u.delta) - B4*x_dot.delta_dot*sin(u.delta) - (A4_dot*sin(u.delta) + A4*x_dot.delta_dot*cos(u.delta));

    % Rear (unsteered) wheels - 
    % wc_2_ax = A2_dot;
    % wc_2_ay = B2_dot;
    % wc_3_ax = A3_dot;
    % wc_3_ay = B3_dot;


% ZeroLongSlipFront
omega_1 = wc_1_vx / P.r_wheel;
omega_4 = wc_4_vx / P.r_wheel;
omega_dot_1 = wc_1_ax / P.r_wheel;
omega_dot_4 = wc_4_ax / P.r_wheel;

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
lambda_1 = (omega_1*P.r_wheel - wc_1_vx) ./ (wc_1_vx + epsVal);
lambda_2 = (omega_2*P.r_wheel - wc_2_vx) ./ (wc_2_vx + epsVal);
lambda_3 = (omega_3*P.r_wheel - wc_3_vx) ./ (wc_3_vx + epsVal);
lambda_4 = (omega_4*P.r_wheel - wc_4_vx) ./ (wc_4_vx + epsVal);

% LateralSlipAngle
alpha_1 = atan(wc_1_vy/wc_1_vx);
alpha_2 = atan(wc_2_vy/wc_2_vx);
alpha_3 = atan(wc_3_vy/wc_3_vx);
alpha_4 = atan(wc_4_vy/wc_4_vx);


Slip = struct();

Slip.long = [lambda_1, lambda_2, lambda_3, lambda_4];
Slip.lat = [alpha_1, alpha_2, alpha_3, alpha_4];


end