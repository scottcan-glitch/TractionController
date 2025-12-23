function u_ref = compute_u_ref_nlp2(xNext_func, x_ref_val, u_guess, u_min, u_max, iterations)
% computes steady-state input for the state x_ref_val.
% iterations is # of forward state calculations the ss input is optimized
% over

    import casadi.*
    
    % delta_dot = 0 for ss
    u_min(4,1) = 0;
    u_max(4,1) = 0;
    
    Nx = numel(x_ref_val);           % state dimension
    Nu = numel(u_guess);              % input dimension

    u = MX.sym('u', Nu);
    xref = MX.sym('xref', Nx);

    
    % steady-state equation error
    x0 = x_ref_val;
    for i=1:iterations
        xNext_val = xNext_func(x0, u);
        x0 = xNext_val;  % Update state for next iteration
    end
    % err = xNext_func(xref, u) - xref;
    err = xNext_val - xref;

    % small regularization for uniqueness & conditioning
    lambda = 1e-6;
    obj = err.'*err + lambda*(u.'*u);

    nlp = struct('x', u, 'f', obj, 'p', xref);
    opts.ipopt.print_level = 0;
    opts.ipopt.tol = 1e-8;

    solver_u = nlpsol('solver_u', 'ipopt', nlp, opts);

    sol = solver_u('x0', u_guess, 'p', x_ref_val, 'lbx', u_min, 'ubx', u_max);
    u_ref = full(sol.x);

end
