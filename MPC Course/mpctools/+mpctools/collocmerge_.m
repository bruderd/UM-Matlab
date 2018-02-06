function vars = collocmerge_(vars, funcargs)
% var = collocmerge_(vars, funcargs)
%
% Merges the 'x' and 'xc' fields of vars. Also merges 'z', and 'zc' if present.
for v = ['x', 'z']
    if isfield(vars, v)
        vars = collocmerge(vars, funcargs, v);
    end
end

end%function

function vars = collocmerge(vars, funcargs, name)
    % name should be a string 'x' or 'z' to indicate which variable is being
    % used.
    namec = [name, 'c'];
    
    if ~all(isfield(vars, {name, namec}))
        error('Variables %s and %s are required for collocation!', ...
              name, namec);
    elseif any(structfun(@(a) ismember(namec, a), funcargs))
        error('%s cannot appear in funcargs (replace with %s)!', ...
              namec, name);
    end
    
    v = vars.(name);
    vc = vars.(namec);
    Nt = length(v) - 1;
    if size(vc, 2) ~= Nt
        error('%s has the wrong number of time points!', namec);
    end
    Nc = size(vc, 1);
    Nv = size(v{1}, 1);

    % Duplicate the endpoints in v and absv.
    vars = rmfield(vars, {name, namec});
    vars.(name) = dupep(v, vc, Nc, Nt);
    absvars = {['abs', name], ['abs', namec]};
    if all(isfield(vars, absvars))
        absv = vars.(absvars{1});
        absvc = vars.(absvars{2});
        vars = rmfield(vars, absvars);
        vars.(absvars{1}) = dupep(absv, absvc, Nc, Nt);
    end
end%function

function X = dupep(x, xc, Nc, Nt)
    % Duplicates the x variables as endpoints of xc.
    X = cell(Nc + 2, Nt);
    for k = 1:Nt
        X(:,k) = [x(k); xc(:,k); x(k + 1)];
    end
end%function

