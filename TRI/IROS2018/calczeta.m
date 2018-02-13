function zeta = calczeta( x, P, params)
%calczeta: Calculates zeta for the parallel combo of FREEs

% Convert end effector to FREE displacement
q = x2q(x, params);

% Calculate the maximum Z = [F,M]' for each FREE
Z = calcZ(q, P, params);

% Convert Z to end effector coordinates
zeta = Z2zeta(Z, params);

end

