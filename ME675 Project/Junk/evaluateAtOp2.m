function f_op = evaluateAtOp2(f, symbols, x_op, u_op)
% evaluateAtOp  Evaluate symbolic model f at a given operating point
%
%   f_op = evaluateAtOp(f, x_op)
%   f_op = evaluateAtOp(f, x_op, u_op)
%
%   f     : symbolic expression / vector / matrix
%   x_op  : numeric or symbolic 5x1 (or 1x5) vector [v_x; v_y; r; omega_2; omega_3]
%   u_op  : (optional) numeric or symbolic 4x1 (or 1x4) vector [T_m2; T_m3; delta; delta_dot]

    % Create canonical symbolic state and input variables (no dependence on workspace)
    x = sym(symbols.State);  % 1x5 sym
    u = sym(symbols.Input.vector, symbols.InputDeriv.delta_dot);    % 1x4 sym

    % Validate inputs
    if nargin < 2
        error('evaluateAtOp requires at least f and x_op.');
    end

    % allow row or column x_op/u_op; convert to column then to row to match x/u
    x_op = sym(x_op);
    x_op = reshape(x_op, 1, []);   % 1x5
    if numel(x_op) ~= numel(x)
        error('x_op must be a 5-element vector (v_x, v_y, r, omega_2, omega_3).');
    end

    if nargin < 3
        do_u = false;
    else
        do_u = true;
        u_op = sym(u_op);
        u_op = reshape(u_op, 1, []);   % 1x4
        if numel(u_op) ~= numel(u)
            error('u_op must be a 4-element vector (T_m2, T_m3, delta, delta_dot).');
        end
    end

    % If f empty -> return immediately
    if isempty(f)
        f_op = f;
        return;
    end

    % Determine which canonical symbols appear in f
    f_syms = symvar(f);    % may return cell or sym array
    if isa(f_syms, 'sym')
        f_sym_names = arrayfun(@char, f_syms, 'UniformOutput', false);
    else
        f_sym_names = cellstr(f_syms);
    end
    x_names = arrayfun(@char, x, 'UniformOutput', false);
    u_names = arrayfun(@char, u, 'UniformOutput', false);

    % Find common canonical state symbols in f (preserves x order)
    [common_x, ia_x, ~] = intersect(x_names, f_sym_names, 'stable');

    % Find common input symbols in f
    [common_u, ia_u, ~] = intersect(u_names, f_sym_names, 'stable');

    if isempty(common_x) && (~do_u || isempty(common_u))
        % Nothing to substitute
        warning('No canonical x or u symbols found in f. Returning f unchanged.');
        f_op = f;
        return;
    end

    % Build substitution lists in the canonical order
    subs_vars = sym([]); subs_vals = sym([]);
    if ~isempty(common_x)
        subs_vars = [subs_vars, x(ia_x)];
        subs_vals = [subs_vals, x_op(ia_x)];
    end
    if do_u && ~isempty(common_u)
        subs_vars = [subs_vars, u(ia_u)];
        subs_vals = [subs_vals, u_op(ia_u)];
    end

    % Perform substitution
    f_op = subs(f, subs_vars, subs_vals);

    % Simplify (optional)
    f_op = simplify(f_op);

    % Warn if still symbolic variables remain
    rem = symvar(f_op);
    if ~isempty(rem)
        if isa(rem, 'sym')
            rem = arrayfun(@char, rem, 'UniformOutput', false);
        else
            rem = cellstr(rem);
        end
        warning('After substitution f_op still contains symbolic variables: %s', strjoin(rem, ', '));
    end
end
