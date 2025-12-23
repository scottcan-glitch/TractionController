function [v_dot_x_out, v_dot_y_out, r_dot_out] = makeBodyFrameDynamics2(State, StateDeriv, Input, P, WheelFrameForces)
% Fast, robust body-frame dynamics solver.
% Tries direct analytic solution if possible; otherwise uses lightweight
% symbolic coefficient extraction or numeric Newton fallback.

% Extract states
v_x = State.v_x;
v_y = State.v_y;
r   = State.r;

% Unknowns (these may be sym or numeric placeholders)
v_dot_x = StateDeriv.v_dot_x;
v_dot_y = StateDeriv.v_dot_y;
r_dot   = StateDeriv.r_dot;

% Parameters
m_vehicle = P.m_vehicle;
I_zz      = P.I_zz;
delta     = Input.delta;

% wheel center positions
wc_1_x = P.wc_1_x; wc_1_y = P.wc_1_y;
wc_2_x = P.wc_2_x; wc_2_y = P.wc_2_y;
wc_3_x = P.wc_3_x; wc_3_y = P.wc_3_y;
wc_4_x = P.wc_4_x; wc_4_y = P.wc_4_y;

% aero
C_drag   = P.C_drag;
rho_air  = P.rho_air;
A_frontal = P.A_frontal;
A_side    = P.A_side;

% Wheel forces in body frame
F_x1 = WheelFrameForces.long(1)*cos(delta) - WheelFrameForces.lat(1)*sin(delta);
F_y1 = WheelFrameForces.lat(1)*cos(delta)  + WheelFrameForces.long(1)*sin(delta);
F_x2 = WheelFrameForces.long(2);
F_y2 = WheelFrameForces.lat(2);
F_x3 = WheelFrameForces.long(3);
F_y3 = WheelFrameForces.lat(3);
F_x4 = WheelFrameForces.long(4)*cos(delta) - WheelFrameForces.lat(4)*sin(delta);
F_y4 = WheelFrameForces.lat(4)*cos(delta)  + WheelFrameForces.long(4)*sin(delta);

F_drag_x = 0.5*C_drag*rho_air*v_x^2*A_frontal;
F_drag_y = 0.5*C_drag*rho_air*v_y^2*A_side;

F_x = F_x1 + F_x2 + F_x3 + F_x4 - F_drag_x;
F_y = F_y1 + F_y2 + F_y3 + F_y4 - F_drag_y;

M_z = wc_1_x * F_y1 - wc_1_y * F_x1 + ...
      wc_2_x * F_y2 - wc_2_y * F_x2 + ...
      wc_3_x * F_y3 - wc_3_y * F_x3 + ...
      wc_4_x * F_y4 - wc_4_y * F_x4;

% Quick check: do the forces depend on the unknown accelerations?
% If not, solution is immediate.
try
    vars_in_F = symvar(F_x) ; % if symbolic, symvar returns symbolic variables in F_x
catch
    vars_in_F = {};
end

dependsOnAccel = true;
if isa(v_dot_x,'sym') || isa(v_dot_y,'sym') || isa(r_dot,'sym')
    % symbolic variables exist. Check if F_x,F_y,M_z contain the unknown symbols.
    dependsOnAccel = any(ismember(vars_in_F, [v_dot_x, v_dot_y, r_dot]));
else
    % numeric path: assume forces are evaluated numerically and do not contain accelerations
    dependsOnAccel = false;
end

% DIRECT analytic solution when F_x,F_y,M_z do NOT depend on accelerations:
if ~dependsOnAccel
    % eqs: v_dot_x - F_x/m - v_y*r == 0  => v_dot_x = F_x/m + v_y*r
    v_dot_x_out = (F_x)/m_vehicle + v_y*r;
    v_dot_y_out = (F_y)/m_vehicle - v_x*r;
    r_dot_out   = (M_z)/I_zz;
    return
end

% If we reach here, F_x/F_y/M_z depend on accelerations (nonlinear/implicit).
% Try to extract linear coefficients quickly (Jacobian) but avoid heavy 'simplify' calls.
% Form residuals as expressions: res = lhs - rhs
res1 = v_dot_x - (F_x/m_vehicle) - v_y*r;
res2 = v_dot_y - (F_y/m_vehicle) + v_x*r;
res3 = r_dot   - (M_z/I_zz);
res = [res1; res2; res3];

% Ensure unknown placeholders are symbolic for jacobian extraction
if ~isa(v_dot_x,'sym') || ~isa(v_dot_y,'sym') || ~isa(r_dot,'sym')
    syms vxd_sym vyd_sym rd_sym real
    res_sym = subs(res, [v_dot_x, v_dot_y, r_dot], [vxd_sym, vyd_sym, rd_sym]);
    z_sym = [vxd_sym; vyd_sym; rd_sym];
else
    res_sym = res;
    z_sym = [v_dot_x; v_dot_y; r_dot];
end

% Compute Jacobian J and RHS B fast (avoid simplify)
J = jacobian(res_sym, z_sym);    % symbolic matrix of partial derivatives
C = subs(res_sym, z_sym, zeros(size(z_sym)));  % constants part
B = -C;

% Try simple symbolic linear solve (J\B). This is fast if J is small and not too complicated.
z_sol = [];
try
    % avoid heavy checks; attempt direct backslash
    z_sol = J \ B;
catch
    z_sol = [];
end

% If symbolic linear solve failed, fallback to numeric Newton solve.
if isempty(z_sol)
    % Numeric fallback: must substitute numeric values for any remaining symbolic parameters.
    % Collect all symbolic parameters needed and try to evaluate numerically.
    % The user should call this function with numeric State/Inputs for simulation.
    try
        % Build numeric function handles for residuals and jacobian
        vars = symvar([res_sym; J; B]);
        % Remove unknown placeholders from vars list; we'll pass unknowns as vector x
        vars = setdiff(vars, z_sym, 'stable');
        % Create matlabFunction for numeric evaluation: F(x,params...)
        Ffun = matlabFunction(res_sym, 'Vars', {z_sym, vars});
        Jfun = matlabFunction(J, 'Vars', {z_sym, vars});
    catch ME
        error('Failed to create numeric function handles for residual/Jacobian: %s', ME.message);
    end

    % Build numeric parameter list from current State/P/Input/P/WheelFrameForces
    % NOTE: We assume here the inputs/states/params are numeric scalars or arrays.
    % Create values vector matching "vars"
    vals = cell(1,numel(vars));
    for i=1:numel(vars)
        name = char(vars(i));
        % Try to extract a value from the workspace of this function's inputs
        % Common symbols will match variable names we have (v_x,v_y,r,omega etc.)
        if isfield(State, name)
            vals{i} = State.(name);
        elseif isfield(StateDeriv, name)
            vals{i} = StateDeriv.(name);
        elseif isfield(Input, name)
            vals{i} = Input.(name);
        elseif isfield(P, name)
            vals{i} = P.(name);
        elseif isfield(WheelFrameForces, name)
            vals{i} = WheelFrameForces.(name);
        else
            % Try evaluating variable from current function workspace (risky)
            try
                vals{i} = eval(name);
            catch
                vals{i} = 0; % fallback; user should provide numeric params
            end
        end
    end

    % initial guess: zeros
    x0 = zeros(3,1);
    opts = optimoptions('fsolve','Display','off','Jacobian','on','FunctionTolerance',1e-8);
    funF = @(x) double(Ffun(x, vals{:}));
    funJ = @(x) double(Jfun(x, vals{:}));

    % Use fsolve with Jacobian approximation (or provide Jacobian)
    try
        xsol = fsolve(@(x) funF(x), x0, opts);
    catch ME
        % final fallback: try vpasolve (slower) or error
        error('Numeric solve via fsolve failed: %s. Consider substituting numeric parameters earlier.', ME.message);
    end

    v_dot_x_out = xsol(1);
    v_dot_y_out = xsol(2);
    r_dot_out   = xsol(3);
    return
end

% z_sol is symbolic (may contain piecewise). Try to return numeric values when possible.
try
    znum = double(z_sol);
    v_dot_x_out = znum(1);
    v_dot_y_out = znum(2);
    r_dot_out   = znum(3);
catch
    % return symbolic expressions
    v_dot_x_out = z_sol(1);
    v_dot_y_out = z_sol(2);
    r_dot_out   = z_sol(3);
end

end
