function testPoints = setTestPoints(testParams, params)
%setTestPoints: Create matrix containing the desired test points
%   Each row is a single test point. Columns are: (1)s, (2)w, (3)P1, (4)P2,
%   (5)P3, ...

Pmax = params.Pmax;
Psteps = testParams.Psteps; % how finely to break up Pmax
stest = testParams.stest;  
wtest = testParams.wtest;

alpha = (0:Psteps)/Psteps;  % fractions of Pmax to take measurements at

Ptest = permn(alpha,3) * diag(Pmax);    % Calulate all test pressure inputs

% enforce minimum pressure, slightly more than 0 psi
Ptest = max(Ptest, ones(size(Ptest))*0.5);

for i = 1:length(stest)
    dim = size(Ptest,1);
    qtest(1:dim , :) = ones(dim,1) * [stest(i), wtest(i)];
    testPoints(dim*(i-1)+1:dim*i , :) = [qtest, Ptest];
end

