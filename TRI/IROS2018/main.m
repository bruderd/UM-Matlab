% main.m

% Set the values of all of the system parameters
Gama = deg2rad([20, -20, 85]);
R = 0.01 * ones(1,3);
L = 0.1 * ones(1,3);
d = zeros(3,3);
p = [0,0,0 ; 0,0,0 ; 1,1,1];
Pmax = [100, 100, 1000];
params = setParams(Gama, R, L, d, p, Pmax);

% Set the value of test parameters
testParams = struct;
testParams.Psteps = 4;     % how finely to break up Pmax
testParams.stest = [0, 5, 5, -5, -5];   % mm
testParams.wtest = [0, 20,-20, 20, -20];    % deg

% Calculate the force zonotope
[zntp, vx, vy] = zonotopeFun(params);

% Create csv of the test data points
testPoints = setTestPoints(testParams, params);
csvwrite('testPoints.csv',testPoints);      % exports testPoint to csv file