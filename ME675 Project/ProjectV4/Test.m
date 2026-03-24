


sedan1 = CarWithDefaults("sedan");
model = Model_v4(sedan1);



for i = 1:numel(model.f)
disp('Sym Vars for')
i
symvar(model.f(i))
end

%%

import casadi.*

% x_op = [1; 1; 0; 1/rad2deg(sedan1.geometry.r_wheel); 1/rad2deg(sedan1.geometry.r_wheel)]
% u_op = [0;0;0;0]

% % Anon. function for getting linear A & B in ct
% A_fun = @(x_op,u_op) evaluateAtOp(model.A,x_op,u_op);
% B_fun = @(x_op,u_op) evaluateAtOp(model.B,x_op,u_op);
% F_fun = @(x_op,u_op) evaluateAtOp(model.f,x_op,u_op);
% F_lin = @(x_op,u_op) F_fun(x_op,u_op) + A_fun(x_op,u_op)*(x-x_op) + B_fun(x_op,u_op)*(u-u_op);


% f_fun = matlabFunction(f_sym,'Vars',{x,u});
x_op = model.symbols.State.vector
u_op = [model.symbols.Input.vector, model.symbols.InputDeriv.delta_dot]

% Continuous time A, B, f at operating point
A_fun = matlabFunction(model.A,'Vars',{x_op,u_op})
B_fun = matlabFunction(model.B,'Vars',{x_op,u_op})
f_fun = matlabFunction(model.f,'Vars',{x_op,u_op})

%% 

% State and Input Sizes
nx = 5;
nu = 4;

% Horizon & Sampling
N = 24; %horizon steps
T = 3;  %horizon time(s)
Ts = T/N;   %time of each step

% x_cl = [x0 zeros(2,Ncost)];
% x_next = x0;
% loops = 20;
% U_optimal = zeros(Ncost,loops);
% u_applied = zeros(1,loops);
% 
% % Initial Condition
% x0 = [0; 0; 0; 0; 0];


% CasADi sym variables
    % State Vector
    x_cas = MX.sym('x_cas',1,nx);
    % v_x = x_cas(1);
    % v_y = x_cas(2);
    % r = x_cas(3);
    % omega_2 = x_cas(4);
    % omega_3 = x_cas(5);

    % Input Vector
    u_cas = MX.sym('u_cas',1,nu);
    % T_m2 = u_cas(1);
    % T_m3 = u_cas(2);
    % delta = u_cas(3);
    % delta_dot = u_cas(4);

% Very light CasADi wrap for calling f_fun
    % F_cas finds F_cas_vector (CasADi symbolic linearized f) at x_cas, u_cas
    % (symbolic CasADi state and input vectors)
if ~( isequal(size(x_cas),[1 5]) && isequal(size(u_cas),[1 4]) )
    disp(x_cas);
    disp(u_cas);
    error('CasADi x_cas or u_cas has wrong dimensions!!');
    
end


f_cas_vector = f_fun([x_cas(1),x_cas(2),x_cas(3),x_cas(4),x_cas(5)],[u_cas(1),u_cas(2),u_cas(3),u_cas(4)])
f_cas = Function('f',{x_cas,u_cas},{f_cas_vector},{'x_cas','u_cas'},{'f_cas'})


% Setup time-integration method. Discretizes the system, numerical
% framework to compute xk+1 from xk, uk. Uses more accurate methods than
% simple forward Euler (1st degree Taylor) discretization

% Integrator to discretize system
intg_options = struct();
intg_options.tf = Ts;
intg_options.simplify = true;
intg_options.number_of_finite_elements = 4;

% Pass ODEs to Integrator (can also be DAEs. Wish I knew this earlier. lol)
dae = struct();
dae.x = x_cas;  %states
dae.p = u_cas;  %parameters (fixed during integration horizon, ie Ts, ZOH)
dae.ode = f_cas(x_cas,u_cas);   %expression for rhs

% Build integrator for discretization, runge-kutta
intg = integrator('intg','rk',dae,intg_options);    % make rk integrator
res = intg('x0',x_cas,'p',u_cas);   % make function for evaluating @ x,u
x_next = res.xf;    % Now return only last element on integration horiz.
F_rk = Function('F_rk',{x_cas,u_cas},{x_next},{'x','u'},{'x_next'});  %wrap, F(x_op,u_op) returns x_k+1

% res = intg('x0',[0;1],'p',u)


% % Forward Euler F_eul instead of rk numerical integration method for
% finding x_k+1. pass

F_fe = @(x_op,u_op) F_eul(x_op,u_op,Ts,A_fun,B_fun,f_fun);

%%


% Open-loop simulation for CasADi rk
Ts = 0.125;
SimTime = 10; %(s)
SimLoops = round(SimTime/Ts);

x_init = [10,0,0,10/sedan1.geometry.r_wheel,10/sedan1.geometry.r_wheel];
u_init = [-100,-100,10,10];
x_now = x_init;
u_now = u_init;

data.x = [x_now; zeros(SimLoops,nx)];
data.u = zeros(SimLoops,nu);

Slip = [straightOnlyGetLongSlip(sedan1,x_init(1,1),x_init(1,4),x_init(1,5)); zeros(SimLoops,4)];
for i=1:SimLoops

    % if i>(SimLoops/2)
    %     u_now = [0,0,0,0];
    % end
x_next = F(data.x(i,:),data.u(i,:));



x_now = x_next;
data.x(i+1,:) = full(x_next);
data.u(i,:) = full(u_now);
Slip(i+1,:) = straightOnlyGetLongSlip(sedan1,data.x(i+1,1),data.x(i+1,4),data.x(i+1,5));
end

% % CHATGPT GARBO

% % -------- Plot States (5) --------
t = (0:SimLoops)' * Ts;   % time vector

figure('Name','States','NumberTitle','off');
for k = 1:6
    if k < 6
        subplot(6,1,k)
        plot(t, data.x(:,k),'LineWidth',1.5)
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
for k = 1:4
    subplot(4,1,k)
    plot(t(1:end-1), data.u(:,k),'LineWidth',1.5)
    grid on
    ylabel(['u_' num2str(k)])
    if k == 1
        title('Input Trajectories')
    end
end
xlabel('Time (s)')


% % CHATGPT GARBO





%%


% Open-loop simulation for Forward-Euler

SimTime = 5; %(s)
Ts = 0.125;
SimLoops = round(SimTime/Ts);

x_init = [10,0,0,10/sedan1.geometry.r_wheel,10/sedan1.geometry.r_wheel];
u_init = [0,0,0,0];
clear time
data.x = [x_init; zeros(SimLoops,nx)];
data.u = zeros(SimLoops,nu);

Slip = [straightOnlyGetLongSlip(sedan1,x_init(1,1),x_init(1,4),x_init(1,5)); zeros(SimLoops,4)];
for i=1:SimLoops

time(i) = i*Ts-Ts
x_now = data.x(i,:);
u_now = u_init;

% % MODIFIED THIS TO ADD x_now 11/16/2025
x_next = F_fe(x_now,u_now)'

data.x(i+1,:) = x_next';
data.u(i,:) = u_now;
Slip(i+1,:) = straightOnlyGetLongSlip(sedan1,data.x(i+1,1),data.x(i+1,4),data.x(i+1,5));
end

% % CHATGPT GARBO

% % -------- Plot States (5) --------
t = (0:SimLoops)' * Ts;   % time vector

figure('Name','States','NumberTitle','off');
for k = 1:6
    if k < 6
        subplot(6,1,k)
        plot(t, data.x(:,k),'LineWidth',1.5)
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
for k = 1:4
    subplot(4,1,k)
    plot(t(1:end-1), data.u(:,k),'LineWidth',1.5)
    grid on
    ylabel(['u_' num2str(k)])
    if k == 1
        title('Input Trajectories')
    end
end
xlabel('Time (s)')


% % CHATGPT GARBO




































%% 

% Optimization problem setup - we can use N concatinated x vectors and N
% concatinated u vectors as the decision variables, subject to system
% dynamics

    % Optimization variables
    opti = casadi.Opti();
    x = opti.variable(nx,N+1);
    u = opti.variable(nu,N);
    % Initial conditions
    x0 = opti.parameter(nx,1);
    delta_dot = opti.parameter(1,N);
    % State cost, input cost, terminal cost >0 (PD). P_cost = none for now
    % Q_cost = 
    % R_cost = 


    % cost function
    opti.minimize(x'*Q_cost*x + u'*R_cost*u)

    % constraints
    opti.subject_to(x(:,1) == x0)
    % opti.subject_to(delta_dot(1) == )
    delta = u(4,:);
    

    for k=1:N
    % dynamics
    opti.subject_to(x(:,k+1) == F_cas(x(:,k),u(:,k)))
    % foward euler for delta_dot
    opti.subject_to((delta(k+1) - delta(k))/Ts == delta_dot(k+1))
    end


% for i=0:loops-1
% 
% 
% 
% 
% 
% Act = A_fun(x_op,u_op);
% Bct = B_fun(x_op,u_op);
% [Ad, Bd] = c2d_exact(A,B,Ts);
% 
% 
% 
% 
% % data.x
% end