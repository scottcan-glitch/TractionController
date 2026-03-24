%% Build HIGH-MU model with sedan1 Car.m instance, sets up CasADi integrator for
% % state evolution. xNext() integrates to next state.

sedan1 = CarWithDefaults_3input("sedan");
model = Model_v5(sedan1);
disp('sym vars?:')
symvar(model.f)

matlabFunction(model.f, 'Vars', {model.symbols.State.vector', model.symbols.StateDeriv.vector', model.symbols.Input.vector' }, 'File', 'dynamicsImplicitFcn');

import casadi.*

Nx = 5;
Nu = 3;
x  = MX.sym('x',Nx,1);       % states
xd = MX.sym('xd',Nx,1);      % state derivatives
u  = MX.sym('u',Nu,1);       % inputs (3x1)

f = dynamicsImplicitFcn(x, xd, u);   % implicit model: 0 = f(x,xd,u, ud)

% Docum: https://web.casadi.org/docs/#nonlinear-programming
dae = struct;
dae.x   = x;      % differential states: the state vector is x
dae.z   = xd;     % treat xdot as algebraic: now we can put the deriv of x into our alg equation (dae.alg)
dae.p   = u;      % parameters/inputs: [T_m2, T_m3, delta]^T
dae.ode = xd;     % trivial ODE equation: the deriv of the state vector is xd
dae.alg = f;      % algebraic constraint form: 0 = f(z,x,p)


% % get next discrete step!!
Ts = 0.0005;  % sample time (fwd integration time)
opts = struct;
opts.tf = Ts;
opts.abstol = 1e-6;
opts.reltol = 1e-6;
opts.max_num_steps = 50000;   % give Newton more space

F = integrator('F', 'idas', dae, opts);

% Wrap integrator setup into callable xNext function
x_k = MX.sym('x_k',Nx,1);
u_k = MX.sym('u_k',Nu,1);

% z_0 = MX.sym('z_0',Nx,1);    %solver guess
% resultFull = F('x0', x_k,'p', u_k,'z0', z_0);   %integrate. z0 is guess for solver
% xNextFull = Function('xNextFull', {x_k, u_k, z_0}, {resultFull.xf, resultFull.zf}); %returns final x_k+1 (Ts later) from integration result

z_0 = MX.zeros(Nx,1);
result = F('x0', x_k,'p', u_k,'z0', z_0);           %integrate. z0 is guess for solver. we could provide a better guess..
xNext = Function('xNext', {x_k, u_k},{result.xf});

% generate c code
xNext.generate('xNext_generated.c')


%% Build Split HIGH/LOW-MU model with sedan2 Car.m instance, sets up CasADi integrator for
% % state evolution. xNext() integrates to next state.

sedan2 = CarWithDefaults_3input("sedan");
sedan2.pacejka.long.mu = [0.3; 0.3; 0.4; 0.4];
sedan2.pacejka.lat.mu = [0.3; 0.3; 0.4; 0.4];

model2 = Model_v5(sedan2);
disp('sym vars?:')
symvar(model2.f)

matlabFunction(model2.f, 'Vars', {model2.symbols.State.vector', model2.symbols.StateDeriv.vector', model2.symbols.Input.vector' }, 'File', 'dynamicsImplicitFcn2');

import casadi.*

Nx = 5;
Nu = 3;
x  = MX.sym('x',Nx,1);       % states
xd = MX.sym('xd',Nx,1);      % state derivatives
u  = MX.sym('u',Nu,1);       % inputs (3x1), last element is delta_dot, which is needed for model.f (due to load xfer = f(delta_dot))

f2 = dynamicsImplicitFcn2(x, xd, u);   % implicit model: 0 = f(x,xd,u, ud)

% Docum: https://web.casadi.org/docs/#nonlinear-programming
dae2 = struct;
dae2.x   = x;      % differential states: the state vector is x
dae2.z   = xd;     % treat xdot as algebraic: now we can put the deriv of x into our alg equation (dae.alg)
dae2.p   = u;      % parameters/inputs: [T_m2, T_m3, delta]^T
dae2.ode = xd;     % trivial ODE equation: the deriv of the state vector is xd
dae2.alg = f2;      % algebraic constraint form: 0 = f(z,x,p)


% % SAME OPTS AS sedan1!
% Ts = 0.001;  % sample time (fwd integration time)
% opts = struct;
% opts.tf = Ts;
% opts.abstol = 1e-6;
% opts.reltol = 1e-6;
% opts.max_num_steps = 50000;   % give Newton more space

F2 = integrator('F2', 'idas', dae2, opts);

% Wrap integrator setup into callable xNext function
x_k = MX.sym('x_k',Nx,1);
u_k = MX.sym('u_k',Nu,1);

% z_0 = MX.sym('z_0',Nx,1);    %solver guess
% resultFull = F('x0', x_k,'p', u_k,'z0', z_0);   %integrate. z0 is guess for solver
% xNextFull = Function('xNextFull', {x_k, u_k, z_0}, {resultFull.xf, resultFull.zf}); %returns final x_k+1 (Ts later) from integration result

z_0 = MX.zeros(Nx,1);
% z_0 = [-0.0741783; -5.97957e-20; 1.07204e-19; 9.52433; 9.52433];
result2 = F2('x0', x_k,'p', u_k,'z0', z_0);   %integrate. z0 is guess for solver
xNext2 = Function('xNext2', {x_k, u_k},{result2.xf});

%% Build LOW-MU model with sedan3 Car.m instance, sets up CasADi integrator for
% % state evolution. xNext() integrates to next state.

sedan3 = CarWithDefaults_3input("sedan");
sedan3.pacejka.long.mu = [0.5; 0.5; 0.5; 0.5];
sedan3.pacejka.lat.mu = [0.5; 0.5; 0.5; 0.5];

model3 = Model_v5(sedan3);
disp('sym vars?:')
symvar(model3.f)

matlabFunction(model3.f, 'Vars', {model3.symbols.State.vector', model3.symbols.StateDeriv.vector', model3.symbols.Input.vector' }, 'File', 'dynamicsImplicitFcn3');

import casadi.*

Nx = 5;
Nu = 3;
x  = MX.sym('x',Nx,1);       % states
xd = MX.sym('xd',Nx,1);      % state derivatives
u  = MX.sym('u',Nu,1);       % inputs (3x1), last element is delta_dot, which is needed for model.f (due to load xfer = f(delta_dot))

f3 = dynamicsImplicitFcn2(x, xd, u);   % implicit model: 0 = f(x,xd,u, ud)

% Docum: https://web.casadi.org/docs/#nonlinear-programming
dae3 = struct;
dae3.x   = x;      % differential states: the state vector is x
dae3.z   = xd;     % treat xdot as algebraic: now we can put the deriv of x into our alg equation (dae.alg)
dae3.p   = u;      % parameters/inputs: [T_m2, T_m3, delta]^T
dae3.ode = xd;     % trivial ODE equation: the deriv of the state vector is xd
dae3.alg = f3;      % algebraic constraint form: 0 = f(z,x,p)


% % SAME OPTS AS sedan1!
% Ts = 0.001;  % sample time (fwd integration time)
% opts = struct;
% opts.tf = Ts;
% opts.abstol = 1e-6;
% opts.reltol = 1e-6;
% opts.max_num_steps = 50000;   % give Newton more space

F3 = integrator('F3', 'idas', dae3, opts);

% Wrap integrator setup into callable xNext function
x_k = MX.sym('x_k',Nx,1);
u_k = MX.sym('u_k',Nu,1);

% z_0 = MX.sym('z_0',Nx,1);    %solver guess
% resultFull = F('x0', x_k,'p', u_k,'z0', z_0);   %integrate. z0 is guess for solver
% xNextFull = Function('xNextFull', {x_k, u_k, z_0}, {resultFull.xf, resultFull.zf}); %returns final x_k+1 (Ts later) from integration result

z_0 = MX.zeros(Nx,1);
% z_0 = [-0.0741783; -5.97957e-20; 1.07204e-19; 9.52433; 9.52433];
result3 = F3('x0', x_k,'p', u_k,'z0', z_0);   %integrate. z0 is guess for solver
xNext3 = Function('xNext3', {x_k, u_k},{result3.xf});
%% BUILD MPC QP ~ CasADi
Ns = 2; % # of slack variables
Nh = 10;
% MPC solver 
solveMPC = makeMPCSolver6(xNext, Nh, sedan1);

    % function returns ss input uref for the reference state xref. this is used
    % in the MPC delta cost function. IE, xNext(xref, uref) ~= xref for steady
    % state. Then x_current, x_ref, uref are passed as params to MPC to
    % calculate u*, which is then injected into system
    u_min = sedan1.limits.Input.Lowerbound;
    u_max = sedan1.limits.Input.Upperbound;
    x_min = sedan1.limits.State.Lowerbound;
    x_max = sedan1.limits.State.Upperbound;
    get_u_ref = @(x_ref,u_guess,xNext) compute_u_ref_ss(xNext,x_ref,u_guess,x_min,x_max,u_min,u_max);

% solveMPC2 = makeMPCSolver(xNext2,sedan2);
%% Test

% % Initial Conditions
% x0_val = [10;0;0;10/sedan1.geometry.r_wheel;10/sedan1.geometry.r_wheel];
% % reference tracking
% xref_val = [40;0;0;40/sedan1.geometry.r_wheel;40/sedan1.geometry.r_wheel];
% uref_val = get_u_ref(xref_val,[26.5;26.5;0;0]);  %if the guess is not good, this may break!
% 
% % numeric parameter vector p0 = [x0_val; xref_val; uref_val]
% p0 = [x0_val; xref_val; uref_val];
% w0 = ones(((Nx+Nu)*Nh),1);%initial guess
% 
% %test
% % sol = solver('x0', w0, 'lbx', args.lbw, 'ubx', args.ubw, 'lbg', args.lbg, 'ubg', args.ubg, 'p', p0);
% % w_opt = full(sol.x)
% 
% disp('DONE!!!!')


%% Simulate!!

Nsim = 2000;

% Initial Condition
x0_val = [28;0;0;28/sedan1.geometry.r_wheel;28/sedan1.geometry.r_wheel];
% x0_val = x_k(:,398); %start from loop 450

% Reference track was28
xref_val = [40;0;0;40/sedan1.geometry.r_wheel;40/sedan1.geometry.r_wheel];
[xref_feasible, uref_feasible] = get_u_ref(xref_val,[100;100;0], xNext)
% Compile into parameter vector p0 = [x0_val; xref_val; uref_val]
p0_val = [x0_val; xref_feasible; uref_feasible];

% Warm Start MPC/rollout (initial solver guess)
    % w0_val = repmat([uref_feasible; xref_feasible],Nh,1);
    % Rollout Version, might work?
    % w0_val = NaN((Nx+Nu)*Nh,1)
    w0_val = [];
    size(w0_val)
    x0 = x0_val;
    S_long_guess = zeros(Ns,1);                                         %this is the slack variable for longitudinal slip. should (almost) always be 0
    for i = 1:Nh
        x0 = xNext(x0, uref_feasible);
        w0_val = [w0_val, [uref_feasible; x0;S_long_guess]];
    end
    disp('Heres w0_val info:')
    w0_val
    size(w0_val)
    w0_val = vertcat(w0_val{:})
    size(w0_val)
    % w0_val = ones(((Nx+Nu)*Nh),1);

% Initialize simulation vectors
x_k = [x0_val, NaN(Nx,Nsim-1)];
u_k = NaN(Nu,Nsim-1);
sol = NaN(Nh);
u_k_horizon = NaN(Nu*Nh,Nsim);
x_k_horizon = NaN(Nx*Nh,Nsim);
s_k_horizon = NaN(Ns*Nh,Nsim);
% status = NaN(Nsim);
lambda = [getLongSlip(x_k(:,1),sedan1), NaN(2,Nsim-1)];

% Bit mask for extracting optimal solution horizons states/inputs/slack
% variables
u_first = [ones(Nu,1); zeros(Nx,1); repmat(zeros(Nu+Nx+Ns,1),Nh-1,1)];
x_first = [zeros(Nu,1); ones(Nx,1); repmat(zeros(Nu+Nx+Ns,1),Nh-1,1)];
u_mask = repmat([ones(Nu,1); zeros(Nx,1); zeros(Ns,1)],Nh,1);
x_mask = repmat([zeros(Nu,1); ones(Nx,1); zeros(Ns,1)],Nh,1);
s_mask = repmat([zeros(Nu,1); zeros(Nx,1); ones(Ns,1)],Nh,1);

%%
% Simulate
for k = 1:1200
    fprintf('Loop Number %d!', k);
    
    % Solve MPC QP
    sol = solveMPC(p0_val,w0_val)
    % Convert soln vector to doubles_k
    w_k = full(sol.x);
    % Record optimal input u* @k
    u_k(:,k) = w_k(u_first == 1);
    % Record k horizon (Nu*Nh & Nx*Nh & Ns*Nh numel columns)
    u_k_horizon(:,k) = w_k(u_mask == 1);
    x_k_horizon(:,k) = w_k(x_mask == 1);
    s_k_horizon(:,k) = w_k(s_mask == 1);
    % Evolve state
    if k<(20)
        % HIGH MU PLANT
        x_kp1 = full(xNext(x_k(:,k),u_k(:,k)));
    
    elseif k >= (20) && k<100
        % Update MPC with the road condition. we assume all states observed, OCP uses the correct slip
        if k == 20
            solveMPC = makeMPCSolver6(xNext2,Nh,sedan2);
        end
        % Split MU PLANT (LOW MU left side, HIGH MU right side)
        x_kp1 = full(xNext2(x_k(:,k),u_k(:,k)));
    elseif k >= (100)
        % Update MPC with the road condition. we assume all states observed, OCP uses the correct slip
        if k == 100
            solveMPC = makeMPCSolver6(xNext3,Nh,sedan3);
        end
        % HIGH MU PLANT
        x_kp1 = full(xNext3(x_k(:,k),u_k(:,k)));
    
    end
    % Guess for next solve
    w0_val = repmat([u_k(:,k); x_kp1;S_long_guess],Nh,1);
    % assign p0 for next MPC solution
    p0_val = [x_kp1; xref_feasible; uref_feasible];

    % Record
    x_k(:,k+1) = x_kp1;
    lambda(:,k+1) = getLongSlip(x_kp1,sedan1);

    %print
    lambda_printout = lambda(:,1:k+1);
    x_k_printout = x_k(:,1:k+1);
    u_k_printout = u_k(:,1:k);
    u_k_horizon_printout = u_k_horizon(:,1:k);
    x_k_horizon_printout = x_k_horizon(:,1:k);
end

%% Plot the results

% time = (Ts * (1:Nsim) - repmat(Ts,Nsim));
% comment below line out to plot everything:
Nsim = 1000;


figure;
subplot(6,1,1);
plot(0:Nsim, x_k(1,1:Nsim+1), 'r-', 'LineWidth', 2);
hold on;
plot(0:Nsim, xref_feasible(1)*ones(1, Nsim+1), 'b--', 'LineWidth', 1.5);
ylim([5,45])
ylabel('Long. Velocity (m/s)');
title('High-Mu Velocity Tracking');
legend('Actual Veloc_x', 'Reference Veloc_x','Location','best');
grid on;


subplot(6,1,2);
plot(0:Nsim, rad2deg(x_k(3,1:Nsim+1)), 'r-', 'LineWidth', 2);
hold on;
plot(0:Nsim, rad2deg(xref_feasible(2)*ones(1, Nsim+1)), 'b--', 'LineWidth', 1.5);
ylim([-5,5])
ylabel('Yaw Rate (deg/s)');
legend('Actual Yaw Rate', 'Reference Yaw Rate','Location','best');
grid on;


% INPUTS
subplot(6,1,3);
stairs(0:Nsim-1, u_k(1,1:Nsim), 'g-', 'LineWidth', 2);
hold on;
plot(0:Nsim, uref_feasible(1)*ones(1, Nsim+1), 'b--', 'LineWidth', 1.5);
ylim([0,200])
xlabel('Time Step');
ylabel('Left Input (Nm)');
title('Control Inputs');
legend('Actual LW Torque', 'Reference LW Torque','Location','best');
grid on;


subplot(6,1,4);
stairs(0:Nsim-1, u_k(2,1:Nsim), 'g-', 'LineWidth', 2);
hold on;
plot(0:Nsim, uref_feasible(2)*ones(1, Nsim+1), 'b--', 'LineWidth', 1.5);
ylim([0,200])
xlabel('Time Step');
ylabel('Right Torque (Nm)');
legend('Actual RW Torque', 'Reference RW Torque','Location','best');
grid on;


subplot(6,1,5);
stairs(0:Nsim-1, rad2deg(u_k(3,1:Nsim)), 'g-', 'LineWidth', 2);
hold on;
plot(0:Nsim, rad2deg(xref_feasible(3)*ones(1, Nsim+1)), 'b--', 'LineWidth', 1.5);
xlabel('Time Step');
ylabel('Steer Angle (deg)');
ylim([-15,15])
legend('Actual Steering Angle', 'Reference Steering Angle','Location','best');
grid on;


subplot(6,1,6);

% Slip
stairs(0:Nsim, 100*lambda(1,1:Nsim+1), 'r--', 'LineWidth', 1.5);
hold on;
stairs(0:Nsim, 100*lambda(2,1:Nsim+1), 'b--', 'LineWidth', 1.5);
ylim([0,24])
xlabel('Time Step (0.5ms)');
ylabel('Long. Slip (%)');
title('Slip');
legend('lambda2', 'lambda3','Location','best');
grid on;

%% Replay Simulation Video

% Example usage assuming you ran your simulation and have the data

% 1. Create a placeholder for the car parameters
% car_params.geometry.a = 1.50; 
% car_params.geometry.b = 1.20;
% car_params.geometry.t = 1.60;
% car_params.geometry.r_wheel = 0.32;
% car_params = sedan1


% 2. Your simulation results (Replace with your actual data)
% X_history: 5xN (Vx, Vy, r, omega2, omega3)
X_history = x_k_printout;
% U_history: 3xN (TL, TR, Delta)
U_history = u_k_printout;
% Ts: Sample time (0.0005s)

create_vehicle_video(X_history, U_history, sedan1, Ts, 'MyAwesomeSim.avi');
