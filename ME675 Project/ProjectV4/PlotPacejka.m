% --- Pacejka Pure Slip Plotting Script ---

% 1. Constants and Parameters
Fz_nominal = 250; % N (Using F_z static value from your car parameters)
eps = 1e-6;        % Small number to prevent division by zero (standard practice)

% Pacejka Parameters (assuming all wheels have the same values)
Blong = 10;
Clong = 1.9;
Elong = 0.97;
Blat  = 8;
Clat  = 1.3;
Elat  = 1.0;

% Coefficients of friction (mu) to plot
mu_values = [0.9, 0.5, 0.2];

% 2. Define Slip Vectors (up to 80% or 0.8 rad)
% Longitudinal Slip (lambda): 0 to 0.8 (80%)
lambda = linspace(-0.8, 0.8, 500); 
% Lateral Slip (alpha): -0.1 to 0.1 radians (~ -5.7 to 5.7 degrees, common working range)
alpha = linspace(-0.1, 0.1, 500); 

% 3. Initialize Figure and Plot Longitudinal Force
figure('Name', 'Pacejka Pure Slip Forces');

% --- Longitudinal Pure Slip (F_long_0 vs lambda) ---
subplot(2, 3, 1:3); % Span across the top row
hold on;
title('Longitudinal Pure Slip Force (F_{x} vs \lambda)');
xlabel('Longitudinal Slip, \lambda');
ylabel('Force (N)');
grid on;

for k = 1:length(mu_values)
    mu_long = mu_values(k);
    
    % The Pure Longitudinal Formula (F_long_0)
    Y = Clong * atan(Blong * lambda - Elong * (Blong * lambda - atan(Blong * lambda)));
    F_long_0 = mu_long * Fz_nominal * sin(Y);
    
    plot(lambda, F_long_0, 'LineWidth', 2, 'DisplayName', ['\mu_{long} = ', num2str(mu_long)]);
end
legend('show', 'Location', 'southeast');
hold off;

% --- Lateral Pure Slip (F_lat_0 vs alpha) ---
subplot(2, 3, 4:6); % Span across the bottom row
hold on;
title('Lateral Pure Slip Force (F_{y} vs \alpha)');
xlabel('Side Slip Angle, \alpha (rad)');
ylabel('Force (N)');
grid on;

for k = 1:length(mu_values)
    mu_lat = mu_values(k);
    
    % The Pure Lateral Formula (F_lat_0)
    Y = Clat * atan(Blat * alpha - Elat * (Blat * alpha - atan(Blat * alpha)));
    F_lat_0 = mu_lat * Fz_nominal * sin(Y);
    
    plot(alpha, F_lat_0, 'LineWidth', 2, 'DisplayName', ['\mu_{lat} = ', num2str(mu_lat)]);
end
legend('show', 'Location', 'southeast');
hold off;
%