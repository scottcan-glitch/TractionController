%% Open loop sim
lambda = [getLongSlip(x_k(:,1),sedan1), NaN(2,Nsim-1)];

Nsim = 20000;
% Initial Condition
x_0 = [10;0;0;10.5/sedan1.geometry.r_wheel;10.5/sedan1.geometry.r_wheel];
x_k = [x_0, NaN(Nx,Nh)];

for k = 1:Nsim
    fprintf('Loop Number %d!', k);

    % Record optimal input u* @k
    u_k(:,k) = [120; 120; deg2rad(0)];

    % Evolve state
    x_kp1 = full(xNext(x_k(:,k),u_k(:,k)));

    % Record
    x_k(:,k+1) = x_kp1
    lambda(:,k+1) = getLongSlip(x_kp1,sedan1);
end

%% Plot the results

% time = (Ts * (1:Nsim) - repmat(Ts,Nsim));

% Fwd velocity
figure;
subplot(6,1,1);
plot(0:Nsim, x_k(1,:), 'r-', 'LineWidth', 2);
hold on;
plot(0:Nsim, xref_feasible(1)*ones(1, Nsim+1), 'b--', 'LineWidth', 1.5);
ylim([5,45])
xlabel('Time Step');
ylabel('Velocity');
title('Velocity Tracking');
legend('Actual Veloc_x', 'Reference Veloc_x');
grid on;

% Lat Velocity
subplot(6,1,2);
plot(0:Nsim, x_k(2,:), 'r-', 'LineWidth', 2);
hold on;
plot(0:Nsim, x_k(3,:), 'b-', 'LineWidth', 2);
xlabel('Time Step');
ylabel('Velocity');
legend('Lateral Velocity', 'Yaw Velocity');
grid on;


% INPUTS
subplot(6,1,3);
stairs(0:Nsim-1, u_k(1,:), 'g-', 'LineWidth', 2);
hold on;
plot(0:Nsim, uref_feasible(1)*ones(1, Nsim+1), 'b--', 'LineWidth', 1.5);
xlabel('Time Step');
ylabel('Left Torque Input');
title('Control Inputs');
legend('Actual LW Torque', 'Reference LW Torque');
grid on;

subplot(6,1,4);
stairs(0:Nsim-1, u_k(2,:), 'g-', 'LineWidth', 2);
hold on;
plot(0:Nsim, uref_feasible(2)*ones(1, Nsim+1), 'b--', 'LineWidth', 1.5);
xlabel('Time Step');
ylabel('Right Torque Input');
legend('Actual RW Torque', 'Reference RW Torque');
grid on;

subplot(6,1,5);
stairs(0:Nsim-1, u_k(3,:), 'g-', 'LineWidth', 2);
hold on;
plot(0:Nsim, xref_feasible(3)*ones(1, Nsim+1), 'b--', 'LineWidth', 1.5);
xlabel('Time Step');
ylabel('Steer Angle Input');
legend('Actual Steering Angle', 'Reference Steering Angle');
grid on;

% Longitudinal Slip
subplot(6,1,6);
plot(0:Nsim, lambda(1,:), 'r--', 'LineWidth', 1.5);
hold on;
plot(0:Nsim, lambda(2,:), 'b--', 'LineWidth', 1.5);
% ylim([0,45])
xlabel('Time Step');
ylabel('Longitudinal Slip');
legend('lambda2', 'lambda3');
grid on;