function solveMPC = makeMPCSolver2(xNext,Nh, sedan)

import casadi.*

Nx = 5;
Nu = 3;
% Makes callable CasADi MPCsolver() which solves NLP based on the dynamics xNext.


    
    
    % Cost matrices (deviation from x_ref, u_ref)
    % Q = diag([60 10 1000 1 1]);   %[fwd_veloc_cost, lat_veloc_cost, yaw_rate_cost, LW_speed_cost, RW_speed_cost]
    % R = diag([1, 1, 20]);  %[LW_torque_cost, RW_torque_cost, steering_angle_cost]
    % P = diag([1 1 1e-3 1 1]);   %[fwd_veloc_cost, lat_veloc_cost, yaw_rate_cost, LW_speed_cost, RW_speed_cost]
    
    Q = diag([50, 400, 1000000, 0.012, 0.012]);
    R = diag([0.007, 0.007, 250000]);
    P = diag([0.2, 800, 2000000, 0.024, 0.024])
    
    % Cost functions
    function l = stageCost(x, xref, u, uref, Q, R)
        l = (x - xref)' * Q * (x - xref) + (u - uref)' * R * (u - uref);
    end
    
    function Vf = termCost(x, xref, P)
        Vf = (x - xref)' * P * (x - xref);
    end
    
    % anonymize
    incCost = @(x, xref, u, uref) stageCost(x,xref,u,uref,Q,R);
    endCost = @(x, xref) termCost(x,xref,P);
    
    % Nh = 10; %horizon
    
    % limit column vectors
    InputLB = sedan.limits.Input.Lowerbound;
    InputUB = sedan.limits.Input.Upperbound;
    StateLB = sedan.limits.State.Lowerbound;
    StateUB = sedan.limits.State.Upperbound;
    
    % syms
    Xk = MX.sym('Xk',Nx,1);  % initial state for solve
    xref = MX.sym('xref',Nx,1);  % reference tracking state
    uref = MX.sym('uref',Nu,1);  % ss input for reference track state
    
    % initialize
    args = struct;
    args.lbw = [];
    args.ubw = [];
    args.lbg = [];
    args.ubg = [];
    args.w0 = []; % init w guess
    nlp = struct;
    nlp.g = {};   % constraints
    nlp.p = vertcat(Xk{:},xref{:},uref{:}); % optimization problem parameters
    nlp.w = {};   %opt_variable
    nlp.J = 0;
    
    % Formulate the NLP
    for k=0:Nh-1
        % New NLP variable for the control input
        Uk = MX.sym(['U_' num2str(k)],Nu,1);
        disp('loop numba!:')
        disp(k)
    
        nlp.w = [nlp.w{:}; Uk] %concat like this for MX vectors. then vertcat(w{:})
        args.lbw = [args.lbw; InputLB]; %torque ref must be under limit, & reference steer angle must be met exactly in a later constraint
        args.ubw = [args.ubw; InputUB];
    
        args.w0 = [args.w0; zeros(Nu,1)] %solver initial guess
    
        % Integrate till the end of the interval
        Xk_end = xNext(Xk, Uk);
        nlp.J = nlp.J + incCost(Xk_end,xref,Uk,uref);
       
    
        % New NLP variable for state at end of interval
        Xk = MX.sym(['X_' num2str(k+1)], Nx,1);
        nlp.w = [nlp.w{:}; Xk];
        args.lbw = [args.lbw; StateLB];
        args.ubw = [args.ubw; StateUB];
        args.w0 = [args.w0; zeros(Nx,1)] %solver initial guess
    
        % Add equality constraint
        nlp.g = [nlp.g, {Xk_end-Xk}];
        args.lbg = [args.lbg; zeros(Nx,1)];
        args.ubg = [args.ubg; zeros(Nx,1)];
    
        % % Enforce: delta = nlp.p{3} = steeringAngleRef
        % these 3 lines were NOT here previously, still testing
        % nlp.g = [nlp.g, { Uk(3) - uref(3) }];
        % args.lbg = [args.lbg; 0];
        % args.ubg = [args.ubg; 0];
    
        % Grab states for ease
        Vx = Xk_end{1};
        r = Xk_end{3};
        omega2 = Xk_end{4};
        omega3 = Xk_end{5};
        r_wheel = sedan.geometry.r_wheel;
        tw = sedan.geometry.t; %trackwidth
    
    
        % Add longitudinal tire slip constraint (<30%)
        wc_2_vx = Vx - (tw/2)*r;
        wc_3_vx = Vx + (tw/2)*r;
        slip_omega2 = (omega2*r_wheel - wc_2_vx) / wc_2_vx;
        slip_omega3 = (omega3*r_wheel - wc_3_vx) / wc_3_vx;
        nlp.g = [nlp.g, {slip_omega2}, {slip_omega3}];
        args.lbg = [args.lbg; -0.30; -0.30];
        args.ubg = [args.ubg; 0.30; 0.30];
    
        % Add lateral tire slip constraint (<20%)
        % wc_2_vy = Vy - dist*r;
        % wc_3_vy = Vy - dist*r;
        % alpha_2 = atan(wc_2_vy/wc_2_vx);
        % alpha_3 = atan(wc_3_vy/wc_3_vx);
        % nlp.g = [nlp.g, {alpha_2}, {alpha_3}];
        % args.lbg = [args.lbg; -0.20; -0.20];
        % args.ubg = [args.ubg; 0.20; 0.20];
    end
    
    % Add terminal cost
    nlp.J = nlp.J + endCost(Xk_end,xref);
    
    nlp_prob = struct('f', nlp.J,'x',vertcat(nlp.w{:}),'g',vertcat(nlp.g{:}),'p',nlp.p);
    
    
    opts2 = struct;
    opts2.ipopt.max_iter = 200;
    opts2.ipopt.print_level = 1; %0,3
    opts2.ipopt.acceptable_tol = 1e-6;
    opts2.ipopt.acceptable_obj_change_tol = 1e-5;
    
    % build solver
    solver = nlpsol('solver', 'ipopt', nlp_prob, opts2);
    
    % wrapper makes easier to use. Only takes params (init x, ref x & u)
    solveMPC = @(p0_initial, w0_guess) solver('x0', w0_guess, 'lbx', args.lbw, 'ubx', args.ubw, 'lbg', args.lbg, 'ubg', args.ubg, 'p', p0_initial);



end