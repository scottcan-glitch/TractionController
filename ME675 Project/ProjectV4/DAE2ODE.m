function f = DAE2ODE(StateDeriv, v_dot_x, v_dot_y, r_dot, omega_dot_2, omega_dot_3)

% Get residual. Second vector term is coupled (DAE), we want explicit ODES
eqs = StateDeriv.vector - [v_dot_x, v_dot_y, r_dot, omega_dot_2, omega_dot_3] == 0;

% Linear solve with Asoln*f = Bsoln
[Asoln,Bsoln] = equationsToMatrix(eqs,StateDeriv.vector);

f = Asoln\Bsoln;

% output is a vector of expressions: 
    % f(1) = expression for v_dot_x
    % f(2) = expression for v_dot_y
    % f(3) = expression for r_dot
    % f(4) = expression for omega_dot_2
    % f(5) = expression for omega_dot_2

end