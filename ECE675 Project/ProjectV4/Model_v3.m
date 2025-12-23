% % ---------------------------------------------------------------------%%
% % The purpose of this script is to develop a nonlinear vehicle dynamics
% % model with 5 states & 3 inputs. The symbolic variables which are intended
% % to be used in the development of a model predictive controller is
% % f_StateFunction, A_nonlinear & B_nonlinear where:
% % 
% %     d(StateVector)/dt = f_StateFunction(StateVector, InputVector)
% %     A_nonlinear = jacobian(f_StateFunction,StateVector)
% %     B_nonlinear = jacobian(f_StateFunction,InputVector)
% % 
% % If A_nonlinear is evaluated with a set of VehicleParameters at the
% % current state and with the current input, the foward euler difference
% % method can be used to approximate system dynamics
% % 
% % VehicleDynamics is the parameter set for the model
% % ---------------------------------------------------------------------%%

function [f_StateFunction, StateVector, InputVector] = Model_v3()

%Steering Angle, Slip Angles, wheel radius
syms delta r_wheel
%Body Frame Velocities
syms v_cg v_x v_y r
%Wheel Center Velocities (in the wheel frame)
syms wc_1_vx wc_1_vy wc_2_vx wc_2_vy wc_3_vx wc_3_vy wc_4_vx wc_4_vy
%Wheel Center Acceleration (in the wheel frame)
syms wc_1_ax wc_1_ay wc_2_ax wc_2_ay wc_3_ax wc_3_ay wc_4_ax wc_4_ay

wc_1_v = [wc_1_vx wc_1_vy]; wc_2_v = [wc_2_vx wc_2_vy];
wc_3_v = [wc_3_vx wc_3_vy]; wc_4_v = [wc_4_vx wc_4_vy];
%Wheel rotational velocity
syms omega_1 omega_2 omega_3 omega_4
%Wheel rotational acceleration
syms omega_dot_1 omega_dot_2 omega_dot_3 omega_dot_4
%Wheel Center locations in Body Frame
syms wc_1 wc_1_x wc_1_y wc_2 wc_2_x wc_2_y wc_3 wc_3_x wc_3_y wc_4 wc_4_x wc_4_y
%Location of COM from wheel contact patches, trackwidth
syms a b h t

BodyFrameWheelLocations = [
wc_1_x == b; wc_1_y == t/2;
wc_2_x == -a; wc_2_y == t/2;
wc_3_x == -a; wc_3_y == -t/2;
wc_4_x == b; wc_4_y == -t/2;
]

WheelFrameWheelVelocity = [
wc_1_vx == (v_x - (wc_1_y*r))*cos(delta) + (v_y + (wc_1_x*r))*sin(delta)
wc_1_vy == (v_y + (wc_1_x*r))*cos(delta) - (v_x - (wc_1_y*r))*sin(delta)
wc_2_vx == v_x - (wc_2_y*r)
wc_2_vy == v_y + (wc_2_x*r)
wc_3_vx == v_x - (wc_3_y*r)
wc_3_vy == v_y + (wc_3_x*r)
wc_4_vx == (v_x - (wc_4_y*r))*cos(delta) + (v_y + (wc_4_x*r))*sin(delta)
wc_4_vy == (v_y + (wc_4_x*r))*cos(delta) - (v_x - (wc_4_y*r))*sin(delta)
]

%Body frame accelerations
syms r_dot v_dot_x v_dot_y

% need to take time deriv of WheelFrameWheelVelocity, use these helper
% functions for product/chain rule, include steer rate terms d(delta)/dt
% Ai, Bi are the terms used in wheel center velocity. ex:
% wc_1_vx = A1*cos(delta) + B1*sin(delta)

    A1 = v_x - wc_1_y*r;   B1 = v_y + wc_1_x*r;
    A2 = v_x - wc_2_y*r;   B2 = v_y + wc_2_x*r;
    A3 = v_x - wc_3_y*r;   B3 = v_y + wc_3_x*r;
    A4 = v_x - wc_4_y*r;   B4 = v_y + wc_4_x*r;
    A1_dot = v_dot_x - wc_1_y*r_dot - wc_1_x*r^2;
    B1_dot = v_dot_y + wc_1_x*r_dot - wc_1_y*r^2;
    A2_dot = v_dot_x - wc_2_y*r_dot - wc_2_x*r^2;
    B2_dot = v_dot_y + wc_2_x*r_dot - wc_2_y*r^2;
    A3_dot = v_dot_x - wc_3_y*r_dot - wc_3_x*r^2;
    B3_dot = v_dot_y + wc_3_x*r_dot - wc_3_y*r^2;
    A4_dot = v_dot_x - wc_4_y*r_dot - wc_4_x*r^2;
    B4_dot = v_dot_y + wc_4_x*r_dot - wc_4_y*r^2;

WheelFrameWheelAcceleration = [
% Front (steered) wheels
wc_1_ax == A1_dot*cos(delta) - A1*sin(delta)*delta_dot + B1_dot*sin(delta) + B1*cos(delta)*delta_dot
wc_1_ay == B1_dot*cos(delta) - B1*delta_dot*sin(delta) - (A1_dot*sin(delta) + A1*delta_dot*cos(delta));

wc_4_ax == A4_dot*cos(delta) - A4*delta_dot*sin(delta) + B4_dot*sin(delta) + B4*delta_dot*cos(delta);
wc_4_ay == B4_dot*cos(delta) - B4*delta_dot*sin(delta) - (A4_dot*sin(delta) + A4*delta_dot*cos(delta));

% Rear (unsteered) wheels
wc_2_ax == A2_dot;
wc_2_ay == B2_dot;
wc_3_ax == A3_dot;
wc_3_ay == B3_dot;
]


%Longitudinal slip
syms lambda_1 lambda_2 lambda_3 lambda_4

ZeroLongSlipFront = [
omega_1 == wc_1_vx / r_wheel
omega_4 == wc_4_vx / r_wheel
omega_dot_1 == wc_1_ax / r_wheel
omega_dot_4 == wc_4_ax / r_wheel
];

ZeroLongSlipFront = [
subs(ZeroLongSlipFront, lhs(WheelFrameWheelVelocity), rhs(WheelFrameWheelVelocity))
]
ZeroLongSlipFront = [
subs(ZeroLongSlipFront, lhs(WheelFrameWheelAcceleration), rhs(WheelFrameWheelAcceleration))
]

%Be careful of divide by zero here, might need max("",1e-6) in den
%also might need abs(wc_4_vx)
LongitudinalSlip = [
lambda_1 == ((omega_1*r_wheel - wc_1_vx)/abs(wc_1_vx))
lambda_2 == ((omega_2*r_wheel - wc_2_vx)/abs(wc_2_vx))
lambda_3 == ((omega_3*r_wheel - wc_3_vx)/abs(wc_3_vx))
lambda_4 == ((omega_4*r_wheel - wc_4_vx)/abs(wc_4_vx))
];

LongitudinalSlip = [
subs(LongitudinalSlip, lhs(ZeroLongSlipFront), rhs(ZeroLongSlipFront))
];
LongitudinalSlip = [
subs(LongitudinalSlip, lhs(WheelFrameWheelVelocity), rhs(WheelFrameWheelVelocity))
]

%Angular slip
syms alpha_1 alpha_2 alpha_3 alpha_4

LateralSlipAngle = [
alpha_1 == atan(wc_1_vy/wc_1_vx)
alpha_2 == atan(wc_2_vy/wc_2_vx)
alpha_3 == atan(wc_3_vy/wc_3_vx)
alpha_4 == atan(wc_4_vy/wc_4_vx)
]

LateralSlipAngle = [
subs(LateralSlipAngle, lhs(WheelFrameWheelVelocity), rhs(WheelFrameWheelVelocity))
]

%Wheel forces
syms F_1_delta F_1_delta_prime F_2_delta F_2_delta_prime F_3_delta F_3_delta_prime F_4_delta F_4_delta_prime
%Wheel forces in body frame
syms F_1_x F_1_y F_2_x F_2_y F_3_x F_3_y F_4_x F_4_y M_1 M_2 M_3 M_4
%Net body frame forces
syms F_x F_y M_z

WheelForceEquations = [
F_1_x == F_1_delta*cos(delta) - F_1_delta_prime*sin(delta)
F_1_y == F_1_delta*sin(delta) + F_1_delta_prime*cos(delta)
F_2_x == F_2_delta
F_2_y == F_2_delta_prime
F_3_x == F_3_delta
F_3_y == F_3_delta_prime
F_4_x == F_4_delta*cos(delta) - F_4_delta_prime*sin(delta)
F_4_y == F_4_delta*sin(delta) + F_4_delta_prime*cos(delta)
]

% Drag Force (assumed x direction only, no moments)
syms F_drag_x F_drag_y C_drag rho_air A_frontal A_side

DragForceEquations = [
F_drag_x == 0.5*C_drag*rho_air*v_x^2*A_frontal;
F_drag_y == 0.5*C_drag*rho_air*v_y^2*A_side
]

% Calculate net forces in the body frame
BodyFrameForceEquations = [
F_x == F_1_x + F_2_x + F_3_x + F_4_x - F_drag_x;
F_y == F_1_y + F_2_y + F_3_y + F_4_y - F_drag_y;
];

BodyFrameForceEquations = [
subs(BodyFrameForceEquations, ...
    [lhs(WheelForceEquations); lhs(DragForceEquations)], ...
    [rhs(WheelForceEquations); rhs(DragForceEquations)])
]

% Calculate the net moment about the z-axis
YawMomentEquation = M_z == wc_1_x * F_1_y - wc_1_y * F_1_x + ... 
wc_2_x * F_2_y - wc_2_y * F_2_x + wc_3_x * F_3_y - wc_3_y * F_3_x + wc_4_x * F_4_y - wc_4_y * F_4_x;

YawMomentEquation = M_z == solve(simplify(subs(YawMomentEquation, ...
    lhs(WheelForceEquations), rhs(WheelForceEquations))), M_z)

%Vehicle Inertias
syms I_zz m_vehicle
% Define the equations of motion for the vehicle
BodyFrameMotionEquations = [
(v_dot_x - v_y*r)*m_vehicle == F_x;
(v_dot_y + v_x*r)*m_vehicle == F_y;
r_dot*I_zz == M_z;
];

BodyFrameMotionEquations = subs(BodyFrameMotionEquations, ...
    [lhs(BodyFrameForceEquations(1)), lhs(BodyFrameForceEquations(2)), lhs(YawMomentEquation)], ...
    [rhs(BodyFrameForceEquations(1)), rhs(BodyFrameForceEquations(2)), rhs(YawMomentEquation)])

%Wheel normal forces using steady state long. and lat. load transfer
%Approximated distribution of roll moment reactive forces between front & rear axle
%lambda_front = ~0.6 sports car understeer, ~0.4 sedan
syms F_z1 F_z2 F_z3 F_z4 lambda_front lambda_rear
syms F_z1_static F_z2_static F_z3_static F_z4_static
syms F_LongitudinalTransfer F_FrontLateralTransfer F_RearLateralTransfer
RollMomentDistribution = [
lambda_front == 0.6
lambda_rear == 1 - lambda_front
]

StaticNormalForce = [
F_z1_static == m_vehicle*9.81*a / (2*(a+b))
F_z2_static == m_vehicle*9.81*b / (2*(a+b))
F_z3_static == m_vehicle*9.81*b / (2*(a+b))
F_z4_static == m_vehicle*9.81*a / (2*(a+b))
]
%Acceleration & left turning loads the right rear (3)
LoadTransfer = [
F_LongitudinalTransfer == m_vehicle*v_dot_x*h/(a+b)
F_FrontLateralTransfer == lambda_front*(m_vehicle*v_dot_y*h/t)
F_RearLateralTransfer == lambda_rear*(m_vehicle*v_dot_y*h/t)
]

LoadTransfer = [
subs(LoadTransfer,lhs(RollMomentDistribution),rhs(RollMomentDistribution))
]
%NormForce = static + long_transfer + lambda*lat_transfer
NormalForce = [
F_z1 == F_z1_static - F_LongitudinalTransfer - F_FrontLateralTransfer;
F_z2 == F_z1_static + F_LongitudinalTransfer - F_RearLateralTransfer;
F_z3 == F_z1_static + F_LongitudinalTransfer + F_RearLateralTransfer;
F_z4 == F_z1_static - F_LongitudinalTransfer + F_FrontLateralTransfer
]
NormalForce = [
subs(NormalForce, lhs(StaticNormalForce), rhs(StaticNormalForce))
]

NormalForce = [
subs(NormalForce, lhs(LoadTransfer), rhs(LoadTransfer))
]

%Power Delivery (total ang. inertias/damping, ang. speed, motor torque)
syms I_1 I_2 I_3 I_4 B_1 B_2 B_3 B_4
syms omega_1 omega_2 omega_3 omega_4 omega_dot_1 omega_dot_2 omega_dot_3 omega_dot_4
syms T_m2 T_m3

PowerTrainEquations = [
-F_1_delta*r_wheel - B_1*omega_1 == I_1*omega_dot_1
T_m2 - F_2_delta*r_wheel - B_2*omega_2 == I_2*omega_dot_2
T_m3 - F_3_delta*r_wheel - B_3*omega_3 == I_3*omega_dot_3
-F_4_delta*r_wheel - B_4*omega_4 == I_4*omega_dot_4
]

PowerTrainEquations = [
subs(PowerTrainEquations,lhs(ZeroLongSlipFront), rhs(ZeroLongSlipFront))
]
PowerTrainEquations = [
subs(PowerTrainEquations,lhs(ZeroLongSlipFront), rhs(ZeroLongSlipFront))
]


%Pacejka MF, friction coeffs., long. & lat.
syms B_MF_long C_MF_long D_MF_long E_MF_long
syms mu_1_long mu_2_long mu_3_long mu_4_long
syms B_MF_lat C_MF_lat D_MF_lat E_MF_lat mu_lat
syms mu_1_lat mu_2_lat mu_3_lat mu_4_lat
Params_Pacejka = [


]
%Pacejka MF-based, wheel 1&4 (front wheels) no longitudinal slip
TractiveEquations = [
F_1_delta == solve(PowerTrainEquations(1), F_1_delta)
F_2_delta == mu_2_long*F_z2*sin(C_MF_long*atan(B_MF_long*lambda_2 - E_MF_long*(B_MF_long*lambda_2 - atan(B_MF_long*lambda_2))))
F_3_delta == mu_3_long*F_z3*sin(C_MF_long*atan(B_MF_long*lambda_3 - E_MF_long*(B_MF_long*lambda_3 - atan(B_MF_long*lambda_3))))
F_4_delta == solve(PowerTrainEquations(4), F_4_delta)

F_1_delta_prime == mu_1_lat*F_z1*sin(C_MF_lat*atan(B_MF_lat*alpha_1 - E_MF_lat*(B_MF_lat*alpha_1 - atan(B_MF_lat*alpha_1))))
F_2_delta_prime == mu_2_lat*F_z2*sin(C_MF_lat*atan(B_MF_lat*alpha_2 - E_MF_lat*(B_MF_lat*alpha_2 - atan(B_MF_lat*alpha_2))))
F_3_delta_prime == mu_3_lat*F_z3*sin(C_MF_lat*atan(B_MF_lat*alpha_3 - E_MF_lat*(B_MF_lat*alpha_3 - atan(B_MF_lat*alpha_3))))
F_4_delta_prime == mu_4_lat*F_z4*sin(C_MF_lat*atan(B_MF_lat*alpha_4 - E_MF_lat*(B_MF_lat*alpha_4 - atan(B_MF_lat*alpha_4))))
]

TractiveEquations = [
subs(TractiveEquations, [lhs(LongitudinalSlip), lhs(LateralSlipAngle), lhs(NormalForce)], ...
[rhs(LongitudinalSlip), rhs(LateralSlipAngle), rhs(NormalForce)])
]

%Compilation into State Space

StateVector = [
v_x
v_y
r
omega_2
omega_3
]

StateVectorDeriv = [
v_dot_x
v_dot_y
r_dot
omega_dot_2
omega_dot_3
]

StateEquations = [
BodyFrameMotionEquations(1)

BodyFrameMotionEquations(2)
BodyFrameMotionEquations(3)
PowerTrainEquations(2)
PowerTrainEquations(3)
]

StateEquations = [
subs(StateEquations,lhs(TractiveEquations),rhs(TractiveEquations))
]

StateEquations = [
v_dot_x == solve(StateEquations(1),v_dot_x)
v_dot_y == solve(StateEquations(2),v_dot_y)
r_dot == solve(StateEquations(3),r_dot)
omega_dot_2 == solve(StateEquations(4),omega_dot_2)
omega_dot_3 == solve(StateEquations(5),omega_dot_3)
]
% The state equation above has coupled derivs of states in the first 3 equations.
% v_dot_x = f(v_dot_y, r_dot), v_dot_y = f(v_dot_x, r_dot), r_dot = f(v_dot_x, v_dot_y)
% Use residual to solve these simultaneously & get each d(statevariable)/dt
% as a function of parameters, inputs & statevariables only (without derivs
% of other state variables involved, which is not ODE form but DAE form,
% where linear algebra is required to simultaneously solve w/ rref etc)

z = lhs(StateEquations(1:3))
Residual = z - rhs(StateEquations(1:3)) == 0
soln = solve(Residual,z)

StateEquations = subs(StateEquations,z,[soln.v_dot_x; soln.v_dot_y; soln.r_dot])
% Now StateEquations is set of 5 coupled ODEs (no longer DAEs)

% x_dot = f_StateFunction(x,u)
f_StateFunction = rhs(StateEquations)

% Input Vector!
InputVector = [
T_m2;
T_m3;
delta
]

% Find A & B matrices for linearization which can be used at each time step.
% Plug in VehicleParameters, StateVector & InputVector to get linear A, B

% A_nonlinear = jacobian(f_StateFunction,StateVector);
% B_nonlinear = jacobian(f_StateFunction,InputVector);


end
% function [] = GetStateFunction()
% 
% 
% end















%Pacejka Parameter List
    %[B_long, B_lat]
    %[C_long, C_lat]
    %[D_long, D_lat]
    %[E_long, E_lat]
%pacejka_ = [10.5, 9.5; 1.65, 1.25; 1.1*F_z, 1.3*F_z; 0.2, -0.5];



% %Vectorized vehicle params
% VehicleParameters = [
% B_1
% B_2
% B_3
% B_4
% B_MF_lat
% B_MF_long
% C_MF_lat
% C_MF_long
% E_MF_lat
% E_MF_long
% I_1
% I_2
% I_3
% I_4
% I_zz
% a
% b
% h
% lambda_front
% m_vehicle
% A_frontal
% A_side
% C_drag
% rho_air
% mu_1_lat
% mu_2_lat
% mu_3_lat
% mu_4_lat
% mu_2_long
% mu_3_long
% r_dot
% r_wheel
% t
% v_dot_x
% v_dot_y
% wc_1_x
% wc_1_y
% wc_2_x
% wc_2_y
% wc_3_x
% wc_3_y
% wc_4_x
% wc_4_y
% wc_4_ax
% ]
% 
% %Solves an expression with input 'Vars'. d/dt(x_bar) in this case
% %combines Ax+Bu (evaluates this with one function of x & u)
% %Takes StateEquations vector of expressions (which symbolically define the
% %time deriv of the state vector), and returns symbolic expression vector
% %evaluated at StateVector, Input, Params
% %d/dt(x_bar) = f_handle(StateVector,Input,VehicleParameters)
% 
% 
% 
% f_handle = matlabFunction(rhs(StateEquations),'Vars',{StateVector,Input,VehicleParameters})
% 
% 
% 
% CasADi
% import casadi.*
% 
% 
% % Symbols/expressions
% x = MX.sym('x');
% y = MX.sym('y');
% z = MX.sym('z');
% f = x^2+100*z^2;
% g = z+(1-x)^2-y;
% 
% nlp = struct;            % NLP declaration
% nlp.x = [x;y;z];         % decision vars
% nlp.f = f;               % objective
% nlp.g = g;               % constraints
% 
% % Create solver instance
% F = nlpsol('F','ipopt',nlp);
% 
% % Solve the problem using a guess
% F('x0',[2.5 3.0 0.75],'ubg',0,'lbg',0)
