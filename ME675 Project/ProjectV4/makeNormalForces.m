function Fz = makeNormalForces(x_dot, car)

% % EXTRACT    
    % Extract car parameters
    m_vehicle = car.dynamics.m_vehicle;
    a = car.geometry.a;
    b = car.geometry.b;
    h = car.geometry.h;
    t = car.geometry.t;
    % Extract RollMomentDistribution
    lambda_front = car.dynamics.lambda_front;
    lambda_rear = 1 - lambda_front;

    % Extract stateDerivs
    v_dot_x = x_dot.v_dot_x;
    v_dot_y = x_dot.v_dot_y;

% % END

% StaticNormalForce


F_static = [m_vehicle*9.81*a/(2*(a+b));
    m_vehicle*9.81*b/(2*(a+b));
    m_vehicle*9.81*b/(2*(a+b));
    m_vehicle*9.81*a/(2*(a+b))];

% LoadTransfer
F_LongitudinalTransfer = m_vehicle*v_dot_x*h/(a+b);
F_FrontLateralTransfer = lambda_front*(m_vehicle*v_dot_y*h/t);
F_RearLateralTransfer = lambda_rear*(m_vehicle*v_dot_y*h/t);

%Acceleration & left turning loads the right rear (3)
%NormForce = static + long_transfer + lambda*lat_transfer
Fz = [
F_static(1) - F_LongitudinalTransfer - F_FrontLateralTransfer;
F_static(2) + F_LongitudinalTransfer - F_RearLateralTransfer;
F_static(3) + F_LongitudinalTransfer + F_RearLateralTransfer;
F_static(4) - F_LongitudinalTransfer + F_FrontLateralTransfer];

end