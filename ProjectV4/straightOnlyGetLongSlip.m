function LongSlip = straightOnlyGetLongSlip(Car, v_x, omega_2,omega_3)

% lambda_1 = (omega_1*r_wheel - wc_1_vx) ./ (wc_1_vx + epsVal);
%  v = r x w
% assumes straight forward motion, such that wc_vx (forward wheel-center velocity)
% is exactly equal to v_x (forward vehicle velocity)

% v_tan_2 = Car.geometry.r_wheel*omega_2;
% v_tan_3 = Car.geometry.r_wheel*omega_3;

epsVal = 1e-6;

LongSlip(1) = 0;
LongSlip(2) = (omega_2*Car.geometry.r_wheel - v_x) ./ (v_x + epsVal);
LongSlip(3) = (omega_3*Car.geometry.r_wheel - v_x) ./ (v_x + epsVal);
LongSlip(4) = 0;

end