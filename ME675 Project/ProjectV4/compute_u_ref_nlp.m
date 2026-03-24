function u_ref = compute_u_ref_nlp(xNext_func, x_ref_val, u_guess, u_min, u_max)

    import casadi.*

    % enforce delta_dot = 0 at steady state
    u_min(4) = 0;
    u_max(4) = 0;

    Nx = numel(x_ref_val);
    Nu = numel(u_guess);

    u     = MX.sym('u',Nu);
    xref  = MX.sym('xref',Nx);

    % steady-state error: xNext(x_ref, u) − x_ref
    xNext_val = xNext_func(xref, u);
    err       = xNext_val - xref;

    % small regularization for uniqueness
    lambda = 0;
    obj = err.'*err + lambda*(u.'*u);

    nlp = struct('x',u,'f',obj,'p',xref);

    opts.ipopt.print_level = 0;
    opts.ipopt.tol         = 1e-8;

    solver_u = nlpsol('solver_u','ipopt',nlp,opts);

    sol   = solver_u('x0',u_guess,...
                     'p',x_ref_val,...
                     'lbx',u_min,...
                     'ubx',u_max);

    u_ref = full(sol.x);

end
