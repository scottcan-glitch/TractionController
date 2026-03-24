function [x_ss_val, u_ss_val] = compute_u_ref_ss(xNext_func, x_ref, u_ref_guess, x_min, x_max, u_min, u_max)

    import casadi.*

    Nx = numel(x_ref);
    Nu = numel(u_ref_guess);

    x_ss = MX.sym('x_ss', Nx);
    u_ss = MX.sym('u_ss', Nu);

    % steady-state constraint: x_next(x_ss, u_ss) = x_ss
    x_next = xNext_func(x_ss, u_ss);
    g = x_next - x_ss;

    % objective: hold x near reference and regularize inputs
    w_x = 1e3;       % penalize deviation from desired state
    w_u = 1e-6;      % small input regularization

    obj = w_x*dot(x_ss - x_ref, x_ss - x_ref) + w_u*dot(u_ss,u_ss);

    nlp = struct('x', [x_ss; u_ss], 'f', obj, 'g', g);

    opts.ipopt.print_level = 0;
    opts.ipopt.tol = 1e-10;

    solver = nlpsol('ss_solver', 'ipopt', nlp, opts);

    x0 = [x_ref; u_ref_guess];
    lbx = [x_min; u_min];
    ubx = [x_max; u_max];

    sol = solver('x0', x0, 'lbx', lbx, 'ubx', ubx, 'lbg', 0, 'ubg', 0);

    x_opt = full(sol.x);
    x_ss_val = x_opt(1:Nx);
    u_ss_val = x_opt(Nx+1:end);

end
