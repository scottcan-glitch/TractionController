function x_next = F_eul(x_op,u_op,Ts,A_fun,B_fun,f_fun)

Ac_lin = A_fun(x_op,u_op);
Bc_linn = B_fun(x_op,u_op);
Bc_lin = [Bc_linn,zeros(5,1)];   %delta_dot is small.
fc_lin = f_fun(x_op,u_op);
I = eye(5);

% change column vectors to vertical
x_op = x_op'
u_op = u_op'
% 
x_now = x_op
u_now = u_op

    % Forward Euler step with affine correction term
    x_next = (I + Ts*Ac_lin)*x_now + Ts*Bc_lin*u_now + Ts*(fc_lin-Ac_lin*x_op-Bc_lin*u_op)
end
