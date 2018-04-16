function Psysid = calcPressureSysidPoints(params)
%setTestPoints: Create matrix containing the desired test points
%   Each row is a single test point. Columns are: (1)P1, (2)P2, (3)P3, ...

num = params.num;   % number of actuators in parallel configuration
pmin = params.pmin;
pmax = params.pmax;
psteps = params.psteps; % how finely to break up Pmax

alpha = (0:psteps)/psteps;  % fractions of Pmax to take measurements at

Ptest = permn(alpha,num) * diag(pmax);    % Calulate all test pressure inputs

% enforce minimum pressure, slightly more than 0 psi
Psysid = max(Ptest, ones(size(Ptest,1),1) * pmin);

end



