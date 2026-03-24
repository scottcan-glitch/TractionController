function lambda = getLongSlip(x_k,sedan1)

% Grab states for ease
Vx = x_k(1,1);
% Vy = x_k(2,1);
r = x_k(3,1);
omega2 = x_k(4,1);
omega3 = x_k(5,1);
r_wheel = sedan1.geometry.r_wheel;
tw = sedan1.geometry.t; %trackwidth
% dist = sedan1.geometry.a; %cg to rear axle

% Add longitudinal tire slip constraint (<30%)
wc_2_vx = Vx - (tw/2)*r;
wc_3_vx = Vx + (tw/2)*r;
lambda_omega2 = (omega2*r_wheel - wc_2_vx) / wc_2_vx;
lambda_omega3 = (omega3*r_wheel - wc_3_vx) / wc_3_vx;

lambda = [lambda_omega2; lambda_omega3];
end