function WheelFrameForces = makeWheelFrameForces(Slip, P, Fz, FrontWheels)
% % Pacejka model used to return wheel frame forces based on wheelFrame
% slip values in the rolling direction (long.) & perpendicular (lat.) for
% the rear wheels, and front wheels for lateral only. Front wheel (1,4)
% longitudinal forces come from powertrain equations (inertia & damping),
% since front wheels assumed a no-slip condition

% long. slip, lat. slip
% lambda_1 = lambda_4 = 0, no slip cond.

% Extract from parameters
B_1 = P.B_1;
B_4 = P.B_4;
I_1 = P.I_1;
I_4 = P.I_4;
r_wheel = P.r_wheel;
% Extract from stateDerivs
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

%pacejka params
for i=1:4

    lambda(i) = Slip.long(i);
    alpha(i) = Slip.lat(i);
    
    Blong(i) = P.pacejka.long.B(i);
    Blat(i) = P.pacejka.lat.B(i);
    Clong(i) = P.pacejka.long.C(i);
    Clat(i) = P.pacejka.lat.C(i);
    Elong(i) = P.pacejka.long.E(i);
    Elat(i) = P.pacejka.lat.E(i);
    mulong(i) = P.pacejka.long.mu(i);
    mulat(i) = P.pacejka.lat.mu(i);

end

% Initialize symbolic vector
    F_long = sym(zeros(1,4));
    F_lat = sym(zeros(1,4));
    F_long_max = sym(zeros(1,4));
    F_lat_max = sym(zeros(1,4));
    s = sym(zeros(1,4));

WheelFrameForces = struct();

    for i=1:4
        % TractiveEquations (Pacejka MF-based, wheel 1&4 (front wheels) no long. slip)
        F_long(i) = mulong(i)*Fz(i)*sin(Clong(i)*atan(Blong(i)*lambda(i) - Elong(i)*(Blong(i)*lambda(i) - atan(Blong(i)*lambda(i)))));
        F_lat(i) = mulat(i)*Fz(i)*sin(Clat(i)*atan(Blat(i)*alpha(i) - Elat(i)*(Blat(i)*alpha(i) - atan(Blat(i)*alpha(i)))));
    
        % max force elipse
        F_long_max(i) = Fz(i)*mulong(i);
        F_lat_max(i) = Fz(i)*mulat(i);
    
        % Scaling factors, if saturated
        s(i) = (F_long(i)/F_long_max(i))^2 + (F_lat(i)/F_lat_max(i))^2;
    
        % Enforce friction elipse using piecewise scaling factor
        F_long(i) = piecewise(s(i)<=1,F_long(i),s(i)>1,F_long(i)/(sqrt(s(i))));
        F_lat(i) = piecewise(s(i)<=1,F_lat(i),s(i)>1,F_lat(i)/(sqrt(s(i))));
    
        % Populate output struct
        WheelFrameForces.long(i) = F_long(i);
        WheelFrameForces.lat(i) = F_lat(i);
    end

    % Front wheel long. forces from drivetrain equ. No slip condition also
    % assumes we stay inside the friction elipse (only lat. is concern for
    % the front wheels to stay inside elipse)

    F_1_long = -(I_1*omega_dot_1 + B_1*omega_1)/r_wheel;
    F_4_long = -(I_4*omega_dot_4 + B_4*omega_4)/r_wheel;

    WheelFrameForces.long(1) = F_1_long;
    WheelFrameForces.long(4) = F_4_long;
end

% F_1_long_max = Fz1*mu_1_long;
% F_1_lat_max = Fz1*mu_1_lat;
% F_2_long_max = Fz2*mu_2_long;
% F_2_lat_max = Fz2*mu_2_lat;
% F_3_long_max = Fz3*mu_3_long;
% F_3_lat_max = Fz3*mu_3_lat;
% F_4_long_max = Fz4*mu_4_long;
% F_4_lat_max = Fz4*mu_4_lat;
% 
% % Scaling factors, if saturated
% s1 = (F_1_long/F_1_long_max)^2 + (F_1_lat/F_1_lat_max)^2;
% s2 = (F_2_long/F_2_long_max)^2 + (F_2_lat/F_2_lat_max)^2;
% s3 = (F_3_long/F_3_long_max)^2 + (F_3_lat/F_3_lat_max)^2;
% s4 = (F_4_long/F_4_long_max)^2 + (F_4_lat/F_4_lat_max)^2;
% 
% % Enforce friction elipse using piecewise scaling factor
% F_1_long = piecewise(s1<=1,F_1_long,s1>1,F_1_long/(sqrt(s1)));
% F_1_lat = piecewise(s1<=1,F_1_lat, s1>1,F_1_lat /(sqrt(s1)));
% F_2_long = piecewise(s2<=1,F_2_long,s2>1,F_2_long/(sqrt(s2)));
% F_2_lat = piecewise(s2<=1,F_2_lat, s2>1,F_2_lat /(sqrt(s2)));
% F_2_long = piecewise(s3<=1,F_3_long,s3>1,F_3_long/(sqrt(s3)));
% F_3_lat = piecewise(s3<=1,F_3_lat, s3>1,F_3_lat /(sqrt(s3)));
% F_4_long = piecewise(s4<=1,F_4_long,s4>1,F_4_long/(sqrt(s4)));
% F_4_lat = piecewise(s4<=1,F_4_lat, s4>1,F_4_lat /(sqrt(s4)));
% 
% 
% 
% WheelFrameForces.long = [F_1_long, F_2_long, F_3_long, F_4_long];
% WheelFrameForces.lat = [F_1_lat, F_2_lat, F_3_lat, F_4_lat];
% 
