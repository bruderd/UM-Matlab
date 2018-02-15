% main.m
clear

%% Set parameter values
% Set the values of system parameters
Gama = deg2rad([40, -40, -89]);
R = (10e-3)/2 * ones(1,3); % R = 0.009525/2 * ones(1,3);
L = 0.1 * ones(1,3);
d = zeros(3,3);
p = [0,0,0 ; 0,0,0 ; 1,1,1];
Pmax = [103.421e3, 103.421e3, 103.421e3];   % Pa = 15psi
params = setParams(Gama, R, L, d, p, Pmax);

% Set the value of test parameters
testParams = struct;
testParams.Psteps = 4;     % how finely to break up Pmax
testParams.stest = [0, 5, 0, 5];   % mm
testParams.wtest = [0, 0, 20, 20];    % deg
testParams.TRmax = 103.421e3;   % Pa = 15psi

% % Much smaller set of points for debugging
% testParams.Psteps = 1;     % how finely to break up Pmax
% testParams.stest = [0, 5];   % mm
% testParams.wtest = [0, 20];    % deg
% testParams.TRmax = 103.421e3;   % Pa = 15psi

%% Calculate the force zonotope for each test configuration
for i = 1:length(testParams.stest)
    s = testParams.stest(i);
    w = testParams.wtest(i);
    x = [0, 0, s*10^(-3), 0, 0, deg2rad(w)]';  % end effector state
    [zntp(:,i), vx(:,i), vy(:,i)] = zonotopeFun(x, params);
end

% Plot all of the force zonotopes using subplots
% ...

%% Create csv of the test data points
testPoints = setTestPoints(testParams, params);
% Need to convert these to voltages before writing to csv or tsv
testPoints_out = testPoints * diag([0.033418, 1, 10/testParams.TRmax, 10/testParams.TRmax, 10/testParams.TRmax]);
csvwrite('testPoints.csv',testPoints_out);      % exports testPoint to csv file
dlmwrite('testPoints.txt',testPoints_out, 'delimiter', '\t', 'newline', 'pc');

%% Calculate the forces at each test data point

predForces = zeros(size(testPoints,1), 2);
for i = 1:size(testPoints,1)
    P = testPoints(i, 3:end)';
    x = [0,0, testPoints(i,1)*10^(-3), 0, 0, deg2rad(testPoints(i,2))]';
    zeta = calczeta(x, P, params);
    predForces(i,:) = [zeta(3), zeta(6)];
end

% Append force predictions to testPoints matrix
testPred = [testPoints, predForces];



