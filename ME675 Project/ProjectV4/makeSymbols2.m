function [x, x_dot, u, u_dot] = makeSymbols2()
% %   Creates a struct of symbolic:
% %         States + State vector
% %         StateDerivs + StateDeriv vector
% %         Input + Input vector
% %         InputDeriv + InputDeriv vector
% %   with dot-access:
% %         x.<IndividualState>
% %         x.vector
% %         x_dot.<IndividualStateDeriv>
% %         x_dot.vector
% %         u.<IndividualInput>
% %         u.vector
% %         u_dot.<IndividualInputDeriv>
% %         u_dot.vector

    % State Vector & Inputs
    stateNames = {'v_x','v_y','r','omega_2','omega_3'};
    stateDerivNames = {'v_dot_x','v_dot_y','r_dot', 'omega_dot_2','omega_dot_3'};
    inputNames = {'T_m2','T_m3','delta'};
    inputDerivNames = {'T_dot_m2', 'T_dot_3', 'delta_dot'};

    % Removes duplicates if any
    stateNames = unique(stateNames, 'stable');
    inputNames = unique(inputNames, 'stable');

% % State Variables & Vector Struct
% % x.r, x.v_x, ... for individual access
% % x.vector for full State Vector
    x = struct();
        for i = 1:numel(stateNames)
            name = stateNames{i};
            % create a symbolic variable and store it as a field (real assumption)
            x.(name) = sym(name, 'real');  
        end

    x.vector = vertcat(cell2sym(cellfun(@(n) x.(n), stateNames, 'UniformOutput', false)));  %vector from struct

% % Symbolic deriv of state vector
    x_dot = struct();
        for i = 1:numel(stateDerivNames)
            name = stateDerivNames{i};
            % create a symbolic variable and store it as a field (real assumption)
            x_dot.(name) = sym(name, 'real');  
        end

    x_dot.vector = vertcat(cell2sym(cellfun(@(n) x_dot.(n), stateDerivNames, 'UniformOutput', false)));  %vector from struct

% % Symbolic deriv of input vector
    u = struct();
        for i = 1:numel(inputNames)
            name = inputNames{i};
            % create a symbolic variable and store it as a field (real assumption)
            u.(name) = sym(name, 'real');  
        end
    u.vector = vertcat(cell2sym(cellfun(@(n) u.(n), inputNames, 'UniformOutput', false)));  %vector from struct

% % Input Variables & Vector Struct
% u.T_m2, u.T_m3, x.delta, ... for individual access
% u.vector for full Input Vector
    u_dot = struct();
        for i = 1:numel(inputDerivNames)
            name = inputDerivNames{i};
            % create a symbolic variable and store it as a field (real assumption)
            u_dot.(name) = sym(name, 'real');  
        end
    u_dot.vector = vertcat(cell2sym(cellfun(@(n) u_dot.(n), inputDerivNames, 'UniformOutput', false)));  %vector from struct
    
end

% % % Descriptors Below:
% % %Steering Angle, Slip Angles, wheel radius
% %   delta r_wheel  
% % %Body Frame Velocities
% %   v_cg v_x v_y r  
% % %Wheel Center Velocities (in the wheel frame)
% %   wc_1_vx wc_1_vy wc_2_vx wc_2_vy wc_3_vx wc_3_vy wc_4_vx wc_4_vy  
% % %Wheel Center Acceleration (in the wheel frame)
% %   wc_1_ax wc_1_ay wc_2_ax wc_2_ay wc_3_ax wc_3_ay wc_4_ax wc_4_ay  
% %   omega_1 omega_2 omega_3 omega_4  
% % %Wheel rotational acceleration
% %   omega_dot_1 omega_dot_2 omega_dot_3 omega_dot_4  
% % %Wheel Center locations in Body Frame
% %   wc_1 wc_1_x wc_1_y wc_2 wc_2_x wc_2_y wc_3 wc_3_x wc_3_y wc_4 wc_4_x wc_4_y  
% % %Location of COM from wheel contact patches, trackwidth
% %   a b h t  
% %   r_dot v_dot_x v_dot_y  
% % %Longitudinal slip
% %   lambda_1 lambda_2 lambda_3 lambda_4  
% % %Angular slip
% %   alpha_1 alpha_2 alpha_3 alpha_4  
% % %Wheel forces
% %   F_1_delta F_1_delta_prime F_2_delta F_2_delta_prime F_3_delta F_3_delta_prime F_4_delta F_4_delta_prime  
% % %Wheel forces in body frame
% %   F_1_x F_1_y F_2_x F_2_y F_3_x F_3_y F_4_x F_4_y M_1 M_2 M_3 M_4  
% % %Net body frame forces
% %   F_x F_y M_z  
% % % Drag Force (assumed x direction only, no moments)
% %   F_drag_x F_drag_y C_drag rho_air A_frontal A_side  
% % %Vehicle Inertias
% %   I_zz m_vehicle  
% % %Wheel normal forces using steady state long. and lat. load transfer
% % %Approximated distribution of roll moment reactive forces between front & rear axle
% % %lambda_front = ~0.6 sports car understeer, ~0.4 sedan
% %   F_z1 F_z2 F_z3 F_z4 lambda_front lambda_rear  
% %   F_z1_static F_z2_static F_z3_static F_z4_static  
% %   F_LongitudinalTransfer F_FrontLateralTransfer F_RearLateralTransfer  
% % %Power Delivery (total ang. inertias/damping, ang. speed, motor torque)
% %   I_1 I_2 I_3 I_4 B_1 B_2 B_3 B_4
% %   omega_1 omega_2 omega_3 omega_4 omega_dot_1 omega_dot_2 omega_dot_3 omega_dot_4  
% %   T_m2 T_m3
% % %Pacejka MF, friction coeffs., long. & lat.
% %   B_long C_long D_long E_long  
% %   mu_1_long mu_2_long mu_3_long mu_4_long  
% %   B_lat C_lat D_lat E_lat mu_lat  
% %   mu_1_lat mu_2_lat mu_3_lat mu_4_lat  


