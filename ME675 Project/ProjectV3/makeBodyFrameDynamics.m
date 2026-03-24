function [v_dot_x, v_dot_y, r_dot] = makeBodyFrameDynamics(State, StateDeriv, Input, P, WheelFrameForces)

% % From this model:
% (v_dot_x - v_y*r)*m_vehicle == F_x;
% (v_dot_y + v_x*r)*m_vehicle == F_y;
% r_dot*I_zz == M_z;

% Extract States
v_x = State.v_x;
v_y = State.v_y;
r = State.r;
v_dot_x = StateDeriv.v_dot_x;
v_dot_y = StateDeriv.v_dot_y;
r_dot = StateDeriv.r_dot;

% % Extract Parameters
m_vehicle = P.m_vehicle;
I_zz = P.I_zz;
delta = Input.delta;

% wheel loc for moments
wc_1_x = P.wc_1_x;
wc_1_y = P.wc_1_y;
wc_2_x = P.wc_2_x;
wc_2_y = P.wc_2_y;
wc_3_x = P.wc_3_x;
wc_3_y = P.wc_3_y;
wc_4_x = P.wc_4_x;
wc_4_y = P.wc_4_y;

% aero
C_drag = P.C_drag;
rho_air = P.rho_air;
A_frontal = P.A_frontal;
A_side = P.A_side;

% Helper expressions
F_x1 = WheelFrameForces.long(1)*cos(delta) - WheelFrameForces.lat(1)*sin(delta);
F_y1 = WheelFrameForces.lat(1)*cos(delta) + WheelFrameForces.long(1)*sin(delta);
F_x2 = WheelFrameForces.long(2);
F_y2 = WheelFrameForces.lat(2);
F_x3 = WheelFrameForces.long(3);
F_y3 = WheelFrameForces.lat(3);
F_x4 = WheelFrameForces.long(4)*cos(delta) - WheelFrameForces.lat(4)*sin(delta);
F_y4 = WheelFrameForces.lat(4)*cos(delta) + WheelFrameForces.long(4)*sin(delta);
F_drag_x = 0.5*C_drag*rho_air*v_x^2*A_frontal;
F_drag_y = 0.5*C_drag*rho_air*v_y^2*A_side;

% Sums
F_x = F_x1 + F_x2 + F_x3 + F_x4 - F_drag_x;
F_y = F_y1 + F_y2 + F_y3 + F_y4 - F_drag_y;
M_z = wc_1_x * F_y1 - wc_1_y * F_x1 + wc_2_x * F_y2 - ...
    wc_2_y * F_x2 + wc_3_x * F_y3 - wc_3_y * F_x3 + ...
    wc_4_x * F_y4 - wc_4_y * F_x4;










% Newtons 3DOF, stored for linear solve w/ z
eq1 = v_dot_x - (F_x/m_vehicle) - v_y*r == 0
eq2 = v_dot_y - (F_y/m_vehicle) + v_x*r == 0
eq3 = r_dot - M_z/I_zz == 0

% [A1, B1] = equationsToMatrix(eq1,v_dot_x);
% eq1_ = v_dot_x - A1\B1 == 0
% [A2, B2] = equationsToMatrix(eq2,v_dot_y);
% eq2_ = v_dot_y - A2\B2 == 0
% [A3, B3] = equationsToMatrix(eq3,r_dot);
% eq3_ = r_dot - A3\B3 == 0

% show classes/sizes that cause the concat error
whos eq1 eq2 eq3 eq1_ eq2_ eq3_
% or
disp('sizes:');
disp(size(eq1)); disp(size(eq2)); disp(size(eq3));
% disp(size(eq1_)); disp(size(eq2_)); disp(size(eq3_));

eqs = [eq1;eq2;eq3];
z_sym = [v_dot_x;v_dot_y;r_dot];
% syms sol1 sol2 sol3
% sol1 == solve(eq1,v_dot_x)
% sol2 == solve(eq2,v_dot_y)
% sol3 == solve(eq3,r_dot)
% eq1_ = v_dot_x - solve(eq1,v_dot_x) == 0
% eq2_ = v_dot_y - solve(eq2,v_dot_y) == 0
% eq3_ = r_dot - solve(eq3,r_dot) == 0





% % GPT

% residuals (make sure these are expressions, not == equalities)
res1 = lhs(eq1) - rhs(eq1);
res2 = lhs(eq2) - rhs(eq2);
res3 = lhs(eq3) - rhs(eq3);
res = [res1; res2; res3];

% unknown vector
z_sym = [v_dot_x; v_dot_y; r_dot];

% 1) Jacobian (coefficients of unknowns)
J = simplify(jacobian(res, z_sym));

% 2) constant/right-hand side: evaluate residuals at z = 0 and negate
B = - simplify(subs(res, z_sym, [0;0;0]));

% display
disp('Jacobian J = ');
disp(J);
disp('B = ');
disp(B);

% 3) Check if J is trivial or singular
detJ = simplify(det(J));

if isAlways(detJ == 0)
    warning('Jacobian determinant is identically zero (symbolic). Singular.');
elseif isAlways(detJ ~= 0)
    % Safe to invert symbolically
    z_sol = simplify(J \ B);
else
    % indeterminate: detJ may be piecewise or depend on parameters
    warning('determinant is indeterminate symbolically (depends on parameters). Attempt symbolic solve and fallback to numeric if needed.');
    try
        z_sol = simplify(J \ B);
    catch
        % fall back to numeric test below or numeric solver
        disp('failed, went to catch')
    end
end



% % GPT
[Asoln,Bsoln] = equationsToMatrix(eqs,z_sym)

z_soln = Asoln\Bsoln


% ONE = size(eq1_)
% TWO = size(eq2_)
% THREE = size(eq3_)


% SOLN = equationsToMatrix([eq1_,eq2_,eq3_],)


% % gpt garbo
% z_sym = [v_dot_x; v_dot_y; r_dot];


% else fall back to case-splitting below




% % gpt garbo end
% 
% % z = [v_dot_x; v_dot_y;r_dot]
% % sys = [solve(eq1,z(1)); solve(eq2,z(2)); solve(eq3,z(3))]
% SOLN = solve([eq1_,eq2_,eq3_],[v_dot_x, v_dot_y,r_dot], 'ReturnConditions',true)
% % sys.
% 
% size_eq1 = size(eq1)
% size_eq2 = size(eq2)
% size_eq3 = size(eq3)
% size_sys = size(sys)
% size_z = size(z)

% Above equations are DAEs;
%   v_dot_x = f(v_dot_y, r_dot)
%   v_dot_y = f(v_dot_x, r_dot)
%   r_dot = f(v_dot_x, v_dot_y)
% Solve residual to reduce to 3 explicit ODEs,
% solve as a system of lin eqs, so that:
%   v_dot_x = f(State, Input, P)
%   v_dot_y = f(State, Input, P)
%   r_dot   = f(State, Input, P)
% Then we will return these ODEs for our final nonlinear f


% [A,B] = equationsToMatrix([eq1;eq2;eq3],z); % returns A*z = B (symbolic)
% 
% z_soln = A\B;
% 
% v_dot_x = z_soln(1);
% v_dot_y = z_soln(2);
% r_dot = z_soln(3);

% z_soln = solve(z-sys == 0, z)




v_dot_x = 1;
v_dot_y = 3;
r_dot = 5;
end