


sedan1 = CarWithDefaults("sedan");
model = Model_v4(sedan1);

matlabFunction(model.f, 'Vars', {model.symbols.State.vector', model.symbols.StateDeriv.vector', [model.symbols.Input.vector, model.symbols.InputDeriv.delta_dot]' }, 'File', 'dynamicsImplicitFcn');

import casadi.*

x  = MX.sym('x',5,1);       % states
xd = MX.sym('xd',5,1);      % state derivatives
u  = MX.sym('u',4,1);       % inputs (3x1), last element is delta_dot, which is needed for model.f (due to load xfer = f(delta_dot))

f = dynamicsImplicitFcn(x, xd, u);   % implicit model: 0 = f(x,xd,u, ud)

% Docum: https://web.casadi.org/docs/#nonlinear-programming
dae = struct;
dae.x   = x;      % differential states: the state vector is x
dae.z   = xd;     % treat xdot as algebraic: now we can put the deriv of x into our alg equation (dae.alg)
dae.p   = u;      % parameters/inputs: [T_m2, T_m3, delta, delta_dot]^T
dae.ode = xd;     % trivial ODE equation: the deriv of the state vector is xd
dae.alg = f;      % algebraic constraint form: 0 = f(z,x,p)


% % get next discrete step!!
Ts = 0.001;  % example sample time (fwd integration time)
opts = struct;
opts.tf = Ts;
opts.abstol = 1e-8;
opts.reltol = 1e-8;
opts.max_num_steps = 50000;   % give Newton more space

F = integrator('F', 'idas', dae, opts);

% Wrap integrator setup into callable xNext function
x_k = MX.sym('x_k',5,1);
u_k = MX.sym('u_k',4,1);
z_0 = MX.sym('z_0',5,1);    %solver guess

result = F('x0', x_k, 'z0', z_0,'p', u_k);   %integrate. z0 is guess for solver
xNext = Function('xNext', {x_k, z_0, u_k}, {result.xf}); %returns final x_k+1 (Ts later) from integration result

% x = [vx; vy; r; omega_2; omega_3]

%% 

% Open Loop Simulation
% x_init = [10;0;0;10/sedan1.geometry.r_wheel;10/sedan1.geometry.r_wheel];
% u_init = [150;150;0;0];
x_init = [10;0;0;10/sedan1.geometry.r_wheel;10/sedan1.geometry.r_wheel];
u_init = [100;101;deg2rad(10);0];
% x_init = [150;0;0;150/sedan1.geometry.r_wheel;150/sedan1.geometry.r_wheel];
% u_init = [0;0;0;0];

nx = 5;
nu = 4;
SimLoops= 10e3;     %for Ts = 0.001, SimLoops = 10e3 --> 10s simulation

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

% xNext()
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
        % ylabel(['x_' num2str(k)])
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
% % State labels
subplot(6,1,1)
ylabel('Long. Veloc(m/s')
subplot(6,1,2)
ylabel('Lat. Veloc(m/s')
subplot(6,1,3)
ylabel('Yaw Rate(rad/s')
subplot(6,1,4)
ylabel('LR Wheel Veloc(rad/s')
subplot(6,1,5)
ylabel('RR Wheel Veloc(rad/s')

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