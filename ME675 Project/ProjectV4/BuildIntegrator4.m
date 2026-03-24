% % Build model with sedan1 Car.m instance, sets up CasADi integrator for
% % state evolution. xNext() integrates to next state.

sedan1 = CarWithDefaults("sedan");
sedan1.geometry.h = 0;
disp('here:')
sedan1.geometry.h
model = Model_v4(sedan1);

matlabFunction(model.f, 'Vars', {model.symbols.State.vector', model.symbols.StateDeriv.vector', [model.symbols.Input.vector, model.symbols.InputDeriv.delta_dot]' }, 'File', 'dynamicsImplicitFcn');

import casadi.*

Nx = 5;
Nu = 4;
x  = MX.sym('x',Nx,1);       % states
xd = MX.sym('xd',Nx,1);      % state derivatives
u  = MX.sym('u',Nu,1);       % inputs (3x1), last element is delta_dot, which is needed for model.f (due to load xfer = f(delta_dot))

f = dynamicsImplicitFcn(x, xd, u);   % implicit model: 0 = f(x,xd,u, ud)

% Docum: https://web.casadi.org/docs/#nonlinear-programming
dae = struct;
dae.x   = x;      % differential states: the state vector is x
dae.z   = xd;     % treat xdot as algebraic: now we can put the deriv of x into our alg equation (dae.alg)
dae.p   = u;      % parameters/inputs: [T_m2, T_m3, delta, delta_dot]^T
dae.ode = xd;     % trivial ODE equation: the deriv of the state vector is xd
dae.alg = f;      % algebraic constraint form: 0 = f(z,x,p)


% % get next discrete step!!
Ts = 0.1;  % example sample time (fwd integration time)
opts = struct;
opts.tf = Ts;
opts.abstol = 1e-8;
opts.reltol = 1e-8;
opts.max_num_steps = 50000;   % give Newton more space

F = integrator('F', 'idas', dae, opts);

% Wrap integrator setup into callable xNext function
x_k = MX.sym('x_k',Nx,1);
u_k = MX.sym('u_k',Nu,1);

% z_0 = MX.sym('z_0',Nx,1);    %solver guess
% resultFull = F('x0', x_k,'p', u_k,'z0', z_0);   %integrate. z0 is guess for solver
% xNextFull = Function('xNextFull', {x_k, u_k, z_0}, {resultFull.xf, resultFull.zf}); %returns final x_k+1 (Ts later) from integration result

z_0 = MX.zeros(Nx,1);
% z_0 = [-0.0741783; -5.97957e-20; 1.07204e-19; 9.52433; 9.52433];
result = F('x0', x_k,'p', u_k,'z0', z_0);   %integrate. z0 is guess for solver
xNext = Function('xNext', {x_k, u_k},{result.xf});

%% BUILD MPC QP ~ CasADi


% function returns ss input uref for the reference state xref. this is used
% in the MPC delta cost function. IE, xNext(xref, uref) ~= xref for steady
% state. Then x_current, x_ref, uref are passed as params to MPC to
% calculate u*, which is then injected into system
u_min = sedan1.limits.Input.Lowerbound;
u_max = sedan1.limits.Input.Upperbound;
x_min = sedan1.limits.State.Lowerbound;
x_max = sedan1.limits.State.Upperbound;
get_u_ref = @(x_ref,u_guess) compute_u_ref_ss(xNext,x_ref,u_guess,x_min,x_max,u_min,u_max);

% Cost matrices
Q = diag([100 100 1e-3 1 1]);
R = diag([0.1, 0.1, 10, 10]);
P = diag([10 10 5 1 1]);

% Cost functions
function l = stageCost(x, xref, u, uref, Q, R)
    l = (x - xref)' * Q * (x - xref) + (u - uref)' * R * (u - uref);
end

function Vf = termCost(x, xref, P)
    Vf = (x - xref)' * P * (x - xref);
end
% anonymize
incCost = @(x, xref, u, uref) stageCost(x,xref,u,uref,Q,R);
endCost = @(x, xref) termCost(x,xref,P);

Nh =3; %horizon

% limit column vectors
InputLB = sedan1.limits.Input.Lowerbound;
InputUB = sedan1.limits.Input.Upperbound;
StateLB = sedan1.limits.State.Lowerbound;
StateUB = sedan1.limits.State.Upperbound;

% syms
Xk = MX.sym('Xk',Nx,1);  % initial state for solve
Uprev = MX.sym('Uprev',Nu,1); % previous input
xref = MX.sym('xref',Nx,1);  % reference tracking state
uref = MX.sym('uref',Nu,1);  % ss input for reference track state

% initialize
args = struct;
args.lbw = [];
args.ubw = [];
args.lbg = [];
args.ubg = [];
args.w0 = []; % init w guess
nlp = struct;
nlp.g = {};   % constraints
nlp.p = vertcat(Xk{:},Uprev{:},xref{:},uref{:}); % optimization problem parameters
nlp.w = {};   %opt_variable
nlp.J = 0;


% Formulate the NLP
for k=0:Nh-1
    % New NLP variable for the control input
    Uk = MX.sym(['U_' num2str(k)],Nu,1);
    disp('loop numba!:')
    disp(k)
    nlp.w = [nlp.w{:}; Uk] %concat like this for MX vectors. then vertcat(w{:})
    args.lbw = [args.lbw; InputLB];
    args.ubw = [args.ubw; InputUB];
    args.w0 = [args.w0;  zeros(Nu,1)] %solver initial guess

    % Integrate till the end of the interval
    Xk_end = xNext(Xk, Uk);
    nlp.J = nlp.J + incCost(Xk_end,xref,Uk,uref); %cost with Xk_end & not Xk, right?
   

    % New NLP variable for state at end of interval
    Xk = MX.sym(['X_' num2str(k+1)], Nx,1);
    nlp.w = [nlp.w{:}; Xk];
    args.lbw = [args.lbw; StateLB];
    args.ubw = [args.ubw;  StateUB];
    args.w0 = [args.w0; zeros(Nx,1)] %solver initial guess

    % Add equality constraint
    nlp.g = [nlp.g, {Xk_end-Xk}];
    args.lbg = [args.lbg; zeros(Nx,1)];
    args.ubg = [args.ubg; zeros(Nx,1)];

    % Enforce: delta_dot = d/dt(delta)
    % over prev applied input & horizon
    % (U_k(3)-U_{k-1}(3))/Ts = U_k(4)
    % if k == 0
    %     % first step uses parameter Uprev (in p0)
    %     prev_u3 = Uprev(3);
    % else
    %     % subsequent steps use previous Uk
    %     prev_u3 = U_prev(3);
    % end
    
    % nlp.g = [nlp.g, { (Uk(3) - prev_u3)/Ts - Uk(4) }];
    % args.lbg = [args.lbg; 0];
    % args.ubg = [args.ubg; 0];

    nlp.g = [nlp.g, {Uk(4)}];
    args.lbg = [args.lbg; 0];
    args.ubg = [args.ubg; 0];


    U_prev = Uk; % store current Uk for next iteration

    % Grab states for ease
    Vx = Xk_end{1};
    Vy = Xk_end{2};
    r = Xk_end{3};
    omega2 = Xk_end{4};
    omega3 = Xk_end{5};
    r_wheel = sedan1.geometry.r_wheel;
    tw = sedan1.geometry.t; %trackwidth
    dist = sedan1.geometry.a; %cg to rear axle

    % Add longitudinal tire slip soft constraint (slack var w/high penalty)
    % when 50% slip is exceeded
    Sk = MX.sym(['S_' num2str(k)],2,1);   % slack for two slips
    nlp.w = [nlp.w{:}; Sk];
    args.lbw = [args.lbw; 0; 0];
    args.ubw = [args.ubw; inf; inf];
    args.w0 = [args.w0; 0; 0];

    wc_2_vx = Vx - (tw/2)*r;
    wc_3_vx = Vx + (tw/2)*r;
    slip_omega2 = (omega2*r_wheel - wc_2_vx) / wc_2_vx;
    slip_omega3 = (omega3*r_wheel - wc_3_vx) / wc_3_vx;
    % nlp.g = [nlp.g, {slip_omega2}, {slip_omega3}];
    % args.lbg = [args.lbg; -0.50; -0.50];
    % args.ubg = [args.ubg; 0.50; 0.50];
    nlp.g = [nlp.g, { slip_omega2 - Sk(1) }, { slip_omega3 - Sk(2) }];
    args.lbg = [args.lbg; -0.50; -0.50];
    args.ubg = [args.ubg;  0.50;  0.50];
    slack_pen = 1e3;  % SLACK PENALTY!
    nlp.J = nlp.J + slack_pen*(Sk(1)^2 + Sk(2)^2);


    % Add lateral tire slip constraint (<20%)
    % wc_2_vy = Vy - dist*r;
    % wc_3_vy = Vy - dist*r;
    % alpha_2 = atan(wc_2_vy/wc_2_vx);
    % alpha_3 = atan(wc_3_vy/wc_3_vx);
    % nlp.g = [nlp.g, {alpha_2}, {alpha_3}];
    % args.lbg = [args.lbg; -0.20; -0.20];
    % args.ubg = [args.ubg; 0.20; 0.20];
end

% Add terminal cost
nlp.J = nlp.J + endCost(Xk_end,xref);

nlp_prob = struct('f', nlp.J,'x',vertcat(nlp.w{:}),'g',vertcat(nlp.g{:}),'p',nlp.p);


opts2 = struct;
opts2.ipopt.max_iter = 200;
opts2.ipopt.print_level = 2; %0,3
opts2.ipopt.acceptable_tol = 1e-6;
opts2.ipopt.acceptable_obj_change_tol = 1e-5;

% build solver
solver = nlpsol('solver', 'ipopt', nlp_prob, opts2);

% wrapper makes easier to use. Only takes params (init x, ref x & u)
solveMPC = @(p0_initial, w0_guess) solver('x0', w0_guess, 'lbx', args.lbw, 'ubx', args.ubw, 'lbg', args.lbg, 'ubg', args.ubg, 'p', p0_initial);


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

Nsim = 20;

% Initial Condition
x0_val = [10;0;0;10/sedan1.geometry.r_wheel;10/sedan1.geometry.r_wheel];
% Reference track
xref_val = [20;0;0;20/sedan1.geometry.r_wheel;20/sedan1.geometry.r_wheel];
disp('attempt ss find for u:')
[xref_val, uref_val] = get_u_ref(xref_val,[26.5;26.5;0;0])
% Compile into parameter vector p0 = [x0_val; xref_val; uref_val]
u_prev_val = [0; 0; 0; 0];  %start of sim, no prev inputs
p0_val = [x0_val; u_prev_val; xref_val; uref_val];

% slack variable for long. slip
S0_val = [0; 0];

% Warm Start MPC (initial solver guess)
w0_val = repmat([uref_val; x0_val; S0_val],Nh,1);
% w0_val = ones(((Nx+Nu)*Nh),1);

% Initialize simulation vectors
x_k = [x0_val, NaN(Nx,Nsim-1)];
u_k = NaN(Nu,Nsim-1);
sol = NaN(Nh);
u_k_horizon = NaN(Nu*Nh,Nsim);
x_k_horizon = NaN(Nx*Nh,Nsim);

% Bit mask for extracting optimal solution states/inputs
u_first = [ones(Nu,1); zeros(Nx,1); repmat(zeros(Nu+Nx,1),Nh-1,1)];
x_first = [zeros(Nu,1); ones(Nx,1); repmat(zeros(Nu+Nx,1),Nh-1,1)];
u_mask = repmat([ones(Nu,1); zeros(Nx,1)],Nh,1);
x_mask = repmat([zeros(Nu,1); ones(Nx,1)],Nh,1);


% Simulate
for k = 1:Nsim
    fprintf('Loop Number %d!', k);
    
    % Solve MPC QP
    sol = solveMPC(p0_val,w0_val)

    % Convert soln vector to double
    w_k = full(sol.x)
    % Record optimal input u* @k
    u_k(:,k) = w_k(u_first == 1)
    u_prev = u_k(:,k);
    % Record k horizon (Nu*Nh & Nx*Nh numel columns)
    u_k_horizon(:,k) = w_k(u_mask == 1)
    x_k_horizon(:,k) = w_k(x_mask == 1)

    % Evolve state
    x_kp1 = full(xNext(x_k(:,k),u_k(:,k)))

    % Guess for next solve
    % shift by one stage: drop first (u,x) and append last stage repeated
    w0_val = [ w_k((Nu+Nx+1):end); w_k(end-(Nu+Nx)+1:end) ];
    % w0_val = repmat([u_k(:,k); x_kp1],Nh,1)

    % assign p0 for next MPC solution
    p0_val = [x_kp1; u_prev; xref_val; uref_val]

    % Record
    x_k(:,k+1) = x_kp1
end

%% Plot the results

% time = (Ts * (1:Nsim) - repmat(Ts,Nsim));

figure;
subplot(2,1,1);
plot(0:Nsim, x_k(1,:), 'r-', 'LineWidth', 2);
hold on;
plot(0:Nsim, xref_val(1)*ones(1, Nsim+1), 'b--', 'LineWidth', 1.5);
xlabel('Time Step');
ylabel('Position');
title('Position Tracking');
legend('Actual Position', 'Reference Position');
grid on;

subplot(2,1,2);
stairs(0:Nsim-1, u_k(1,:), 'g-', 'LineWidth', 2);
xlabel('Time Step');
ylabel('Control Input');
title('Control Input Over Time');
grid on;

