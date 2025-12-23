function J = f_cost(U,N_horizon,Q,R,P)
% U is a vertically concatenated vector of control input vectors u.
% Our system is single input so ui is 1x1
% For L input system & N horizon, U has dimension LN x 1. Function
% returns cost for the series of inputs over the horizon (applies u1,
% u2, ..., uN)

J = 0;
    for k=0:N_horizon-1
        u_inc = U(k+1);
        J_inc = sum(u_inc.^2);
        J = J + J_inc;
    end
end