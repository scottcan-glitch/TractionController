%% Minimal edits of your script to use xNextFull + disable calc_ic + warm zf
clear; clc;
sedan1 = CarWithDefaults("sedan");
model = Model_v4(sedan1);

matlabFunction(model.f, 'Vars', {model.symbols.State.vector', model.symbols.StateDeriv.vector', [model.symbols.Input.vector, model.symbols.InputDeriv.delta_dot]' }, 'File', 'dynamicsImplicitFcn');

import casadi.*

Nx = 5;
Nu = 4;
x  = MX.sym('x',Nx,1);       % states
xd = MX.sym('xd',Nx,1);      % state derivatives
u  = MX.sym('u',Nu,1);       % inputs

f = dynamicsImplicitFcn(x, xd, u);   % implicit model: 0 = f(x,xd,u)

dae = struct;
dae.x   = x;      
dae.z   = xd;     % treat xdot as algebraic
dae.p   = u;      
dae.ode = xd;     
dae.alg = f;      

% ------------------ IMPORTANT: disable calc_ic so IDAS won't call IDACalcIC during jacobian/sens calls
Ts = 0.1;
opts = struct;
opts.tf = Ts;
opts.abstol = 1e-8;
opts.reltol = 1e-8;
opts.max_num_steps = 50000;
opts.calc_ic = false;   % <-- Minimal edit #1: prevent IDAS from trying to compute consistent ICs during sensitivity calls

F = integrator('F', 'idas', dae, opts);

% ------------------ Create xNextFull that returns both xf and zf (uncommented/used)
x_k = MX.sym('x_k',Nx,1);
u_k = MX.sym('u_k',Nu,1);
z_0 = MX.sym('z_0',Nx,1);    % solver guess for algebraic variables

resFull = F('x0', x_k, 'p', u_k, 'z0', z_0);
xNextFull = Function('xNextFull', {x_k, u_k, z_0}, {resFull.xf, resFull.zf}, ...
                     {'x','u','z0'}, {'xfinal','zfinal'});

% Build mpctools wrapper that takes x,u and an extra parameter z0 (z0 will be provided to nmpc as a parameter)
mpc = import_mpctools();   % you said you fixed instance issue earlier
xNextFull_mp = mpc.getCasadiFunc(@(xx,uu,zz0) xNextFull(xx,uu,zz0), [Nx,Nu,Nx], {'x','u','z0'}, {'xplus'});

% -------- Costs (unchanged)
Q = diag([10 10 1e-3 1e-3 1e-3]);
R = diag([.1, .1, 1, 1]);
P = diag([10 10 5 5 1]);
xref = zeros(5,1);

l = mpc.getCasadiFunc(@(x,u) stageCost(x,u,xref,Q,R), [Nx, Nu], {'x', 'u'}, {'l'});
Vf = mpc.getCasadiFunc(@(x) termCost(x,xref,P), Nx, {'x'}, {'Vf'});

% -------- MPC sizes & guesses (unchanged except we provide par.z0)
Nt = 2;        % Prediction horizon
Nsim = 10;

N = struct('x', Nx,'u',Nu,'t',Nt);

x0 = [10;0;0;10/sedan1.geometry.r_wheel;10/sedan1.geometry.r_wheel];
u0 = [100;100;0;0];

% initial z0 numeric guess for algebraic states (zeros is OK if reasonable)
z0_init = zeros(Nx,1);

% Build initial rollout guess using the explicit integrator call (we use xNext (now we also have xNextFull available))
X0 = repmat(x0,1,Nt+1);
U0 = repmat(u0,1,Nt);
x_tmp = x0;

for kk=1:Nt
    % Use xNextFull here to be consistent: pass z0_init as guess
    [x_tmp, z_tmp] = xNextFull_mp(x_tmp, u0, z0_init);
    x_tmp = full(x_tmp);
    X0(:,kk+1) = x_tmp;
end
guess = struct('x',X0,'u',U0);

% Create parameter struct for NMPC: par.z0 is expected by the wrapped function (one column per horizon step)
par = struct();
par.z0 = repmat(z0_init, 1, Nt);   % Minimal edit #2: provide z0 parameter as time-varying across horizon

% -------- Build solver ONCE (do not rebuild in loop)
solver = mpc.nmpc('f', xNextFull_mp, 'N', N, 'l', l, 'x0', x0, 'guess', guess, 'par', par, 'Vf', Vf, 'verbosity', '3');

% MUST initialize solver before calling solve()
solver.init();

% -------- Simulation loop (minimal edits)
x = [x0, zeros(Nx,Nsim)];
u = [u0, zeros(Nu,Nsim-1)];

x_k = x0;
u_k = u0;
zf_prev = z0_init;   % carry zf from plant integration to use as warm guess

for k = 1:Nsim
    disp(['step ' num2str(k) ':'])

    % Reset internal solver state (you used this earlier)
    solver.reset();

    % Fix initial state at first horizon index (1-based)
    solver.fixvar('x',1,x_k);

    % Update parameter z0 for this iteration with the latest algebraic guess repeated across horizon
    par.z0 = repmat(zf_prev, 1, Nt);
    % Minimal edit #3: set par into solver via cyclepar (mpctools API)
    solver.cyclepar('z0', par.z0);

    % Solve (no args)
    solver.solve();
    disp(['solver.status: ' solver.status])

    % Check results
    if isequal(solver.status, 'Solve_Succeeded') || isequal(solver.status, 'Solved_To_Acceptable_Level')
        solver.saveguess();   % warm-start for next iter
    else
        warning('Solver failed at time step %d!',k);
        break
    end

    % get first control
    u_k = full(solver.var.u(:,1));
    disp('u_k = '); disp(u_k')

    % Integrate the true plant using xNextFull and supply zf_prev as z0 to get consistent zf
    [x_next_sym, zf_sym] = xNextFull(x_k, u_k, zf_prev);  % returns MX; evaluate below
    x_kp1 = full(x_next_sym);
    zf_prev = full(zf_sym);

    disp('x_k+1 = '); disp(x_kp1')

    % store results
    u(:,k) = u_k;
    x(:,k+1) = x_kp1;
    x_k = x_kp1;

    % DO NOT rebuild solver here — we keep the same solver
end

% (unchanged) plotting or further processing can go here

%% cost functions
function l = stageCost(x, u, xref, Q, R)
error = x - xref;
l = transpose(error) * Q * error + transpose(u) * R * u;
end

function Vf = termCost(x, xref, P)
error = x - xref;
Vf = transpose(error) * P * error;
end
