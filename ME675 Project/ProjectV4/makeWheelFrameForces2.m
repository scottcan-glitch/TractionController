function WheelFrameForces = makeWheelFrameForces2(Slip, car, Fz, FrontWheels)
% % Pacejka model used to return wheel frame forces based on wheelFrame
% slip values in the rolling direction (long.) & perpendicular (lat.) for
% the rear wheels, and front wheels for lateral only. Front wheel (1,4)
% longitudinal forces come from powertrain equations (inertia & damping),
% since front wheels assumed a no-slip condition

% long. slip, lat. slip
% lambda_1 = lambda_4 = 0, no slip cond.

% No divide by zero.
eps = 1e-4;
% % EXTRACT
    % Extract from car parameters
    B_1 = car.dynamics.B_1;
    B_4 = car.dynamics.B_4;
    I_1 = car.dynamics.I_1;
    I_4 = car.dynamics.I_4;
    r_wheel = car.geometry.r_wheel;
    % Extract front wheel speeds, solved from kinematics (no slip)
    omega_dot_1 = FrontWheels.omega_dot_1;
    omega_dot_4 = FrontWheels.omega_dot_4;
    omega_1 = FrontWheels.omega_1;
    omega_4 = FrontWheels.omega_4;

% % Initialize symbolics
    lambda = sym(zeros(1,4));
    alpha = sym(zeros(1,4));
    Blong = sym(zeros(1,4));
    Blat  = sym(zeros(1,4));
    Clong = sym(zeros(1,4));
    Clat  = sym(zeros(1,4));
    Elong = sym(zeros(1,4));
    Elat  = sym(zeros(1,4));
    mulong = sym(zeros(1,4));
    mulat  = sym(zeros(1,4));

    % Extract pacejka params
    for i=1:4
    
        lambda(i) = Slip.long(i);
        alpha(i) = Slip.lat(i);
        
        Blong(i) = car.pacejka.long.B(i);
        Blat(i) = car.pacejka.lat.B(i);
        Clong(i) = car.pacejka.long.C(i);
        Clat(i) = car.pacejka.lat.C(i);
        Elong(i) = car.pacejka.long.E(i);
        Elat(i) = car.pacejka.lat.E(i);
        mulong(i) = car.pacejka.long.mu(i);
        mulat(i) = car.pacejka.lat.mu(i);
    
    end

    % % END

% Initialize symbolic vector
    F_long_0 = sym(zeros(1,4));
    F_lat_0 = sym(zeros(1,4));
    F_long = sym(zeros(1,4));
    F_lat = sym(zeros(1,4));

WheelFrameForces = struct();

% No Slip Front wheels
WheelFrameForces.long(1) = -(I_1*omega_dot_1 + B_1*omega_1)/r_wheel;
WheelFrameForces.long(4) = -(I_4*omega_dot_4 + B_4*omega_4)/r_wheel;

    for i=1:4

        % TractiveEquations (Pacejka MF-based, wheel 1&4 (front wheels) no long. slip)
        if i == 2 || i == 3
            F_long_0(i) = mulong(i)*Fz(i)*sin(Clong(i)*atan(Blong(i)*lambda(i) - Elong(i)*(Blong(i)*lambda(i) - atan(Blong(i)*lambda(i)))));
            F_lat_0(i) = mulat(i)*Fz(i)*sin(Clat(i)*atan(Blat(i)*alpha(i) - Elat(i)*(Blat(i)*alpha(i) - atan(Blat(i)*alpha(i)))));
        
            % Scale for smooth combined friction results
            F_lat(i) = F_lat_0(i)/(eps+sqrt(1+(F_long_0(i)/mulong(i))^2));
            F_long(i) = F_long_0(i)/(eps+sqrt(1+(F_lat_0(i)/mulat(i))^2));
        elseif i == 1 || i == 4
            % F_long_0(i) = mulong(i)*Fz(i)*sin(Clong(i)*atan(Blong(i)*lambda(i) - Elong(i)*(Blong(i)*lambda(i) - atan(Blong(i)*lambda(i)))));
            F_lat_0(i) = mulat(i)*Fz(i)*sin(Clat(i)*atan(Blat(i)*alpha(i) - Elat(i)*(Blat(i)*alpha(i) - atan(Blat(i)*alpha(i)))));
            F_long(i) = WheelFrameForces.long(i);
            F_long_0(i) = F_long(i)*(eps+sqrt(1+(F_lat_0(i)/mulat(i))^2));

            % Scale for smooth combined friction results
            F_lat(i) = F_lat_0(i)/(eps+sqrt(1+(F_long_0(i)/mulong(i))^2));
            % F_long(i) = F_long_0(i)/(sqrt(1+(F_lat_0(i)/mulat(i)^2)));

            % F_long(i) = WheelFrameForces.long(i)
        end

    WheelFrameForces.long(i) = F_long(i);
    WheelFrameForces.lat(i) = F_lat(i);

    end

    % Front wheel long. forces from drivetrain equ. No slip condition also
    % assumes we stay inside the friction elipse (only lat. is concern for
    % the front wheels to stay inside elipse)

end
