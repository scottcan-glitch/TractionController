


sedan1 = CarWithDefaults("sedan");
model = Model_v4(sedan1);



for i = 1:numel(model.f)
disp('Sym Vars for')
i
symvar(model.f(i))
end

% states and stuff
% syms v_x v_y r omega_2 omega_3 real
% syms v_dot_x v_dot_y r_dot omega_2_dot omega_3_dot real
% syms T_m2 T_m3 delta real
% 
% x = [v_x; v_y; r; omega_2; omega_3];
% u = [T_m2; T_m3; delta];
% x_dot = [v_dot_x; v_dot_y; r_dot; omega_2_dot; omega_3_dot];


matlabFunction(model.f, 'Vars', {model.symbols.State.vector', model.symbols.StateDeriv.vector', model.symbols.Input.vector'}, 'File', 'dynamicsImplicitFcn');

% % % % % % % % % % setup in casadi
% import casadi.*
% 
% x = MX.sym('x',1,5);      % states
% xd = MX.sym('xd',1,5);    % state derivatives
% u = MX.sym('u',1,3);      % inputs
% 
% f = dynamicsImplicitFcn(x, xd, u);   % call MATLAB-generated function
% 
% dae = struct;
% dae.x = x;      % differential states
% dae.z = xd;     %algebraic
% dae.p = u;      % inputs treated as parameters
% dae.ode = f;    % implicit dynamics

% % % % % % % % % % % % 
import casadi.*

x  = MX.sym('x',5,1);       % states
xd = MX.sym('xd',5,1);      % state derivatives
u  = MX.sym('u',3,1);       % inputs

f = dynamicsImplicitFcn(x, xd, u);   % implicit model: 0 = f(x,xd,u)

% Docum: https://web.casadi.org/docs/#nonlinear-programming
dae = struct;
dae.x   = x;      % differential states: the state vector is x
dae.z   = xd;     % treat xdot as algebraic: now we can put the deriv of x into our alg equation (dae.alg)
dae.p   = u;      % parameters/inputs:         
dae.ode = xd;     % trivial ODE equation: the deriv of the state vector is xd
dae.alg = f;      % algebraic constraint form: 0 = f(z,x,p)


% % get next discrete step!!
Ts = 0.001;  % example sample time
opts = struct;
opts.tf = Ts;
opts.abstol = 1e-8;
opts.reltol = 1e-8;
opts.max_num_steps = 50000;   % give Newton more space

F = integrator('F', 'idas', dae, opts);

% Wrap integrator setup into callable xNext function
x_k = MX.sym('x_k',5,1);
u_k = MX.sym('u_k',3,1);

result = F('x0', x_k, 'z0', zeros(5,1),'p', u_k);   %integrate
xNext = Function('xNext', {x_k, u_k}, {result.xf}); %returns final x state (Ts later) from integration


% x = [vx; vy; r; omega_2; omega_3]

%% 

% Simulation Setup
x_init = [10;0;0;10/sedan1.geometry.r_wheel;10/sedan1.geometry.r_wheel];
u_init = [150;150;0];
% x_init = [15;15;0;15/sedan1.geometry.r_wheel;15/sedan1.geometry.r_wheel];
% u_init = [10;10;0];

nx = 5;
nu = 3;
SimLoops= 4000;

data.x = zeros(nx,SimLoops+1);
data.u = zeros(nu,SimLoops);

data.x(:,1) = x_init;
data.u = repmat(u_init,1,SimLoops);


Slip = [straightOnlyGetLongSlip(sedan1,x_init(1),x_init(4),x_init(5)); zeros(SimLoops,4)];
% disp('thing')
% size(Slip)

xd_prev = zeros(5,1);

for i=1:SimLoops

x_op = data.x(:,i);
u_op = data.u(:,i);


result = F('x0', data.x(:,i), 'z0', xd_prev, 'p', data.u(:,i));
data.x(:,i+1) = full(result.xf);
xd_prev = full(result.zf);%store for next loop guess

% data.x(:,i+1) = full(xNext(x_op,u_op));
Slip(i+1,:) = straightOnlyGetLongSlip(sedan1,data.x(1,i+1),data.x(4,i+1),data.x(5,i+1));

% % 
end

% % CHATGPT GARBO

% % -------- Plot States (5) --------
t = (0:SimLoops)' * Ts;   % time vector

figure('Name','States','NumberTitle','off');
for k = 1:6
    if k < 6
        subplot(6,1,k)
        plot(t, data.x(k,:),'LineWidth',1.5)
        grid on
        ylabel(['x_' num2str(k)])
    end
    if k == 6
        subplot(6,1,k)
        plot(t,Slip(:,2:3),'LineWidth',1.5)
        title('Longitudinal Slip Ratio')
    end
    if k == 1
        title('State Trajectories')
    end
end
xlabel('Time (s)')

% % -------- Plot Inputs (4) --------
figure('Name','Inputs','NumberTitle','off');
for k = 1:3
    subplot(4,1,k)
    plot(t(1:end-1), data.u(k,:),'LineWidth',1.5)
    grid on
    ylabel(['u_' num2str(k)])
    if k == 1
        title('Input Trajectories')
    end
end
xlabel('Time (s)')


% % CHATGPT GARBO

%%

import casadi.*

x  = MX.sym('x',5,1);
xd = MX.sym('xd',5,1);
u  = MX.sym('u',3,1);

f = dynamicsImplicitFcn(x, xd, u);

dae = struct;
dae.x   = x;
dae.z   = xd;
dae.p   = u;
dae.ode = xd;
dae.alg = f;

Ts = 0.01;
F  = integrator('F','idas',dae,struct('tf',Ts));

% initial condition
x_init = [10;0;0;10/sedan1.geometry.r_wheel;10/sedan1.geometry.r_wheel];
u_init = [100;100;0];

SimLoops = 200;
nx = 5; nu = 3;

data.x = zeros(nx,SimLoops+1);
data.u = repmat(u_init,1,SimLoops);
data.x(:,1) = x_init;

for i = 1:SimLoops
    res = F('x0', data.x(:,i), 'z0', zeros(5,1), 'p', data.u(:,i));
    x_next = full(res.xf);
    data.x(:,i+1) = x_next;
end

% % -------- Plot Inputs (4) --------
figure('Name','Inputs','NumberTitle','off');
for k = 1:3
    subplot(4,1,k)
    plot(t(1:end-1), data.u(k,:),'LineWidth',1.5)
    grid on
    ylabel(['u_' num2str(k)])
    if k == 1
        title('Input Trajectories')
    end
end
xlabel('Time (s)')