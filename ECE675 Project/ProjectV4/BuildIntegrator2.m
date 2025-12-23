
% % Builds model with sedan1 Car.m instance, sets up CasADi integrator for
% % state evolution. xNext() integrates to next state.

sedan1 = CarWithDefaults("sedan");
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

%% BUILD MPC

addpath('C:\Users\123sc\MPCTools\rawlings-group-octave-mpctools-765025c73156')
mpc = import_mpctools();

% Ts must be the same in this script as BuildIntegrator.m, make this
% robust!!

Q = diag([10 10 1e-3 1e-3 1e-3]);
R = diag([.1, .1, 1, 1]);
P = diag([10 10 5 5 1]);
% xref = [10;0;0;10/sedan1.geometry.r_wheel;10/sedan1.geometry.r_wheel];
xref = zeros(5,1);

% x_ref = MX.sym('x_ref',Nx,N+1);
% u_prev = MX.sym('u_prev',Nu);


Nt = 2;        % Prediction horizon
Nsim = 10;    % Sim steps. 10e3 = 10s (for Ts = 0.001s)
% U = MX.sym('U',Nu,Nsim);       % control increments Δu
% X = MX.sym('X',Nx,Nsim+1);     % predicted states


% lb.u = double(sedan1.limits.Input.Lowerbound);
% ub.u = double(sedan1.limits.Input.Upperbound);
% lb.x = double(sedan1.limits.State.Lowerbound);
% ub.x = double(sedan1.limits.State.Upperbound);
N = struct('x', Nx,'u',Nu,'t',Nt);
% N.c = 2;  %collocation
% get all the CasADi funcs

% l = mpc.getCasadiDAE(F)
l = mpc.getCasadiFunc(@(x,u) stageCost(x,u,xref,Q,R), [Nx, Nu], {'x', 'u'}, {'l'});
Vf = mpc.getCasadiFunc(@(x) termCost(x,xref,P), Nx, {'x'}, {'Vf'});
xNext_mp = mpc.getCasadiFunc(@(xx,uu) xNext(xx,uu), [Nx,Nu],{'x','u'},{'xplus'});
% l: stagecost, lb:lowerbound, ub:upperbound, Vf: terminalcost, f:
% stateevolution (discrete next state from x, u) N: size of u x & horizon,
% kwargs: holds/passes a bunch of these args, verbosity: solver tells you
% more info with higher int args




x0 = [10;0;0;10/sedan1.geometry.r_wheel;10/sedan1.geometry.r_wheel];
u0 = [100;100;0;0];
x = [x0, zeros(Nx,Nsim)];
u = [u0, zeros(Nu,Nsim-1)];


% % compute rollout for solver initial guesses
X0 = repmat(x0,1,Nt+1);
U0 = repmat(u0,1,Nt);
x_tmp = x0;
for kk=1:Nt
    x_tmp = full( xNext(x_tmp, u0) );   % or xNext_mp
    X0(:,kk+1) = x_tmp;
end

guess = struct('x',X0,'u',U0)
% solver = mpc.nmpc('f',xNext_mp,'N',N,'l',l,'x0',x0,'lb',lb,'ub',ub,'Vf',Vf,'verbosity','3','guess',guess);

opts = struct('tol', 1e-6,'acceptable_tol',1e-4);
% opts.tol = 1e-6;              % Relax tolerance if too strict
% opts.acceptable_tol_ = 1e-4;   % Acceptable solution tolerance
% opts.max_iter = 500;          % Increase iterations if needed
% opts.linear_solver = 'ma57';  % Or 'mumps' if available
% opts.hessian_approximation = 'limited-memory'; % Faster, avoids full Hessian
% Delta=Ts*5
solver = mpc.nmpc('f',xNext_mp,'N',N,'l',l,'x0',x0,'guess',guess);

% solver.guess('u', repmat(u0,1,Nt));

% initial points
x_k = x0;
u_k = u0;
    for k = 1:Nsim
        disp('step1!')
        
        solver.reset()
        % solver.fixvar('x',1,x_k) %sets initial value for nmpc horizon
        % solver.fixvar('u',1,u_k)
        disp('step2!')
        solver.solve()  %solves nmpc
        solver.status
        disp('step3!')
        if isequal(solver.status, 'Solve_Succeeded')
            solver.saveguess();
            % solver.fixvar('x',1,solver.var.x(:,1));
        elseif isequal(solver.status, 'Solved_To_Acceptable_Level')
            solver.saveguess();
            warning('Solved_To_Acceptable_Level at time step %d!',k)
        else
            warning('Solver failed at time step %d!',k);
            break
        end
        disp('step4!')
        u_k = full(solver.var.u(:,1))
        x_kp1 = full(xNext(x_k,u_k))
        thing = size(x_kp1)
        disp('step5!')
        % solver.fixvar('x',1,x_kp1)

        u(:,k) = u_k;
        x(:,k+1) = x_kp1;
        
        %update guess?
        guess = struct('x',x_kp1,'u',u_k)
        solver = mpc.nmpc('f',xNext_mp,'N',N,'l',l,'x0',x_kp1,'verbosity','0','guess',guess)
    end


% mpc.mpcplot(x, u, t)

function l = stageCost(x, u, xref, Q, R)
error = x - xref;
l = transpose(error) * Q * error + transpose(u) * R * u;
% l = transpose(error) * error + transpose(u) * u;
end

function Vf = termCost(x, xref, P)
error = x - xref;

Vf = transpose(error) * P * error;
% Vf = transpose(error) * error;
end

rmpath('C:\Users\123sc\MPCTools\rawlings-group-octave-mpctools-765025c73156')