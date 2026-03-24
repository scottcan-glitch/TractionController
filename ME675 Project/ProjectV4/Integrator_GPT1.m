%% NMPC example using implicit model + CasADi idas integrator
clear; clc;
import casadi.*;

%% Build model & integrator
sedan1 = CarWithDefaults("sedan");
model = Model_v4(sedan1);

matlabFunction(model.f, 'Vars', {model.symbols.State.vector', model.symbols.StateDeriv.vector', [model.symbols.Input.vector, model.symbols.InputDeriv.delta_dot]' }, 'File', 'dynamicsImplicitFcn');

Nx = 5;
Nu = 4;

x  = MX.sym('x',Nx,1);
xd = MX.sym('xd',Nx,1);
u  = MX.sym('u',Nu,1);

f = dynamicsImplicitFcn(x, xd, u);

dae = struct;
dae.x   = x;
dae.z   = xd;
dae.p   = u;
dae.ode = xd;
dae.alg = f;

Ts = 0.1;
opts = struct;
opts.tf = Ts;
opts.abstol = 1e-8;
opts.reltol = 1e-8;
opts.max_num_steps = 50000;

F = integrator('F', 'idas', dae, opts);

x_k = MX.sym('x_k',Nx,1);
u_k = MX.sym('u_k',Nu,1);
result = F('x0', x_k,'p', u_k,'z0', zeros(Nx,1));
xNext = Function('xNext', {x_k, u_k},{result.xf});

%% --- MPC SETUP ---
addpath('C:\Users\123sc\MPCTools\rawlings-group-octave-mpctools-765025c73156')
mpc = import_mpctools();

Q = diag([10 10 1e-3 1e-3 1e-3]);
R = diag([.1 .1 1 1]);
P = diag([10 10 5 5 1]);

xref = zeros(Nx,1);

Nt = 12;      % prediction horizon
Nsim = 100;   % closed-loop sim steps

l  = mpc.getCasadiFunc(@(x,u) stageCost(x,u,xref,Q,R), [Nx, Nu], {'x','u'}, {'l'});
Vf = mpc.getCasadiFunc(@(x) termCost(x,xref,P), Nx, {'x'}, {'Vf'});
xNext_mp = mpc.getCasadiFunc(@(xx,uu) xNext(xx,uu), [Nx,Nu], {'x','u'}, {'xplus'});

N = struct('x',Nx,'u',Nu,'t',Nt);

%% INITIAL GUESS TRAJECTORIES
x0 = [10;0;0;10/sedan1.geometry.r_wheel;10/sedan1.geometry.r_wheel];
u0 = [100;100;0;0];

X0 = repmat(x0,1,Nt+1);
U0 = repmat(u0,1,Nt);

x_tmp = x0;
for kk = 1:Nt
    x_tmp = full(xNext(x_tmp,u0));
    X0(:,kk+1) = x_tmp;
end

guess = struct('x',X0,'u',U0);

%% BUILD SOLVER ONCE
solver = mpc.nmpc('f',xNext_mp,'N',N,'l',l,'Vf',Vf,'x0',x0,'guess',guess,'verbosity',3);

%% CLOSED LOOP SIMULATION
x_k = x0;
x    = zeros(Nx,Nsim+1);
u    = zeros(Nu,Nsim);
x(:,1) = x_k;

for k = 1:Nsim
    fprintf('\n--- MPC STEP %d ---\n', k)

    solver.fixvar('x',1,x_k);   % set initial state

    solver.solve();
    disp(solver.status);

    if isequal(solver.status,'Solve_Succeeded') || isequal(solver.status,'Solved_To_Acceptable_Level')
        solver.saveguess();
        solver.shift();
    else
        warning('Solver failure at k=%d, resetting guess...',k)
        solver.reset();
        continue;
    end

    u_k = full(solver.var.u(:,1));
    x_kp1 = full(xNext(x_k,u_k));

    u(:,k)   = u_k;
    x(:,k+1) = x_kp1;
    x_k      = x_kp1;
end

%% Plot
t = 0:Ts:Nsim*Ts;
figure; plot(t,x'); title('State trajectory')
figure; stairs(t(1:end-1),u'); title('Input trajectory')

rmpath('C:\Users\123sc\MPCTools\rawlings-group-octave-mpctools-765025c73156')

%% costs
function l = stageCost(x, u, xref, Q, R)
err = x - xref;
l = err'*Q*err + u'*R*u;
end

function Vf = termCost(x, xref, P)
err = x - xref;
Vf = err'*P*err;
end
