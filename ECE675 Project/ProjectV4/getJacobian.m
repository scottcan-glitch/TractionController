function jacob = getJacobian(f, q)
% getJacobian  Return jacobian of f with respect to q.vector
%   f - symbolic expression (or vector)
%   q - struct with field 'vector' that contains symbolic variables
%
%   jacob - symbolic Jacobian (or [] if variables missing)

    % default return if we exit early
    jacob = sym([]);

    % validate input
    if ~isstruct(q) || ~isfield(q, 'vector')
        error('Input q must be a struct with a field ''vector''.');
    end

    % ensure q_vector is a symbolic row vector
    q_vector = sym(q.vector);             % convert to symbolic if not already
    q_vector = reshape(q_vector, 1, []);  % force row vector

    % get symbol names as cell arrays of chars for reliable comparison
    f_symobjs = symvar(f);            % may return sym array or cell array
    q_symobjs = symvar(q_vector);

    % convert sym (or char) entries to cell array of character vectors
    if isempty(f_symobjs)
        f_symbols = {};
    elseif isa(f_symobjs, 'sym')
        f_symbols = arrayfun(@char, f_symobjs, 'UniformOutput', false);
    else
        % e.g., symvar may already return a cell array of char
        f_symbols = cellstr(f_symobjs);
    end

    if isempty(q_symobjs)
        q_symbols = {};
    elseif isa(q_symobjs, 'sym')
        q_symbols = arrayfun(@char, q_symobjs, 'UniformOutput', false);
    else
        q_symbols = cellstr(q_symobjs);
    end

    % If q_vector contains no symbolic variables, return jacobian of f wrt empty -> empty
    if isempty(q_symbols)
        warning('q.vector contains no symbolic variables. Returning empty Jacobian.');
        return;
    end

    % check all q_symbols appear in f
    if ~all(ismember(q_symbols, f_symbols))
        warning('Some elements of q.vector do not appear in f. Returning empty Jacobian.');
        return;
    end

    % compute jacobian (works if f is scalar, vector, or matrix)
    jacob = jacobian(f, q_vector);
end
