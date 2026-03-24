addpath('C:\Users\123sc\MPCTools\rawlings-group-octave-mpctools-765025c73156')
mpc = import_mpctools();
import casadi.*

% Ts must be the same in this script as BuildIntegrator.m, make this
% robust!!

Q = diag([10 10 5 5 1]);
R = diag([1 1 0 0]);
P = diag([10 10 5 5 1]);
xref = [12;0;0;12/sedan1.geometry.r_wheel;12/sedan1.geometry.r_wheel];

% x_ref = MX.sym('x_ref',Nx,N+1);
% u_prev = MX.sym('u_prev',Nu);

Nx = 5;
Nu = 4;
Nt = 20;        % Prediction horizon
Nsim = 10e3;    % Sim steps. 10e3 = 10s (for Ts = 0.001s)
% U = MX.sym('U',Nu,Nsim);       % control increments Δu
% X = MX.sym('X',Nx,Nsim+1);     % predicted states


lb.u = sedan1.limits.Input.Lowerbound;
ub.u = sedan1.limits.Input.Upperbound;
lb.x = sedan1.limits.State.Lowerbound;
ub.x = sedan1.limits.State.Upperbound;
N = struct('x', Nx,'u',Nu,'t',Nt);

% get all the CasADi funcs

% l = mpc.getCasadiDAE(F)
l = mpc.getCasadiFunc(@(x,u) stageCost(x,u,xref,Q,R), [Nx, Nu], {'x', 'u'}, {'l'});
Vf = mpc.getCasadiFunc(@(x) termCost(x,xref,P), Nx, {'x'}, {'Vf'});
xNext = mpc.getCasadiFunc(@(x,u) xNext(x,u,z_0), [Nx,Nu],{'x','u'});
% l: stagecost, lb:lowerbound, ub:upperbound, Vf: terminalcost, f:
% stateevolution (discrete next state from x, u) N: size of u x & horizon,
% kwargs: holds/passes a bunch of these args, verbosity: solver tells you
% more info with higher int args
kwargs = struct('l',l,'lb',lb,'ub',ub,'Vf',Vf,'verbosity','3');
solver = mpc.nmpc('f',xNext,'N',N,'**',kwargs);



x0 = [10;0;0;10/sedan1.geometry.r_wheel;10/sedan1.geometry.r_wheel];
u0 = [100;101;deg2rad(10);0];


    for i = 1:Nsim

        
        solver.fixvar('x',0,x0) %sets initial value for nmpc horizon
        solver.solve()  %solves nmpc
        if isequal(solver.status, 'Solve_Succeeded')
            solver.saveguess()
        else
            warning('Solver failed at time step %d!',ii);
            break
        end


    end


% mpc.mpcplot(x, u, t)

function l = stageCost(x, u, xref, Q, R)
error = x - xref;
l = transpose(error) * Q * error + transpose(u) * R * u;
end

function Vf = termCost(x, xref, P)
error = x - xref;

Vf = transpose(error) * P * error;
end

rmpath('C:\Users\123sc\MPCTools\rawlings-group-octave-mpctools-765025c73156')