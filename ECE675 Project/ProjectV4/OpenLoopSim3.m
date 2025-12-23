%% SIM

% Initial Conditions
x_init = [10;0;0;10/sedan1.geometry.r_wheel;10/sedan1.geometry.r_wheel];
u_init = [100;101;deg2rad(10);0];

% Sim time
SimLoops= 10e3;     %for Ts = 0.001, SimLoops = 10e3 --> 10s simulation

% Data and OL input
nx = size(x_init,1);    %5
nu = size(u_init,1);    %4
data.x = zeros(nx,SimLoops+1);
data.u = zeros(nu,SimLoops);
data.x(:,1) = x_init;
data.u = repmat(u_init,1,SimLoops); %open loop constant input
z_guess = zeros(5,1);

% Slip data
Slip = [straightOnlyGetLongSlip(sedan1,x_init(1),x_init(4),x_init(5)); zeros(SimLoops,4)];

for i=1:SimLoops

% Get current point and input
x_k = data.x(:,i);
u_k = data.u(:,i);

% Integrated to next point Ts later
[x_kp1, z_guess] = xNext(x_k,u_k,z_guess);    %returns x_k+1, and start point guess for next integration.
    



% Record state
data.x(:,i+1) = full(x_kp1);

% Record slip
Slip(i+1,:) = straightOnlyGetLongSlip(sedan1,data.x(1,i+1),data.x(4,i+1),data.x(5,i+1));

end

%% PLOT

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