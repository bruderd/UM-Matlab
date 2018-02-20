function [FM1, FM2, FM3, FM4, FM1pred, FM2pred, FM3pred, FM4pred, RMSE] = treatData(MatDataFile, testPoints2, testPoints_out, params, testParams)
% treatData: Reads in measured data, removes elastomer offsets, splits into
% separate arrays for each configuration. Also measures error
% MatDataFile is a string which is the name of the .mat file containing the
% test data

% Load in test data and input points
load(MatDataFile);
FMraw = SyncFT.data;
RESraw = SycRes.data;

% remove bad transition points (where linear actuator or servo weren't set)
i = 1;
while i <= length(RESraw)
    if (RESraw(i,7) == 125 || RESraw(i,7) == 145) && abs((RESraw(i,8)-4.077)*(1/0.033418) - 2.5) > 2
        i = i+1;
    else
        RESraw(i,:) = [];
        FMraw(i,:) = [];
    end
end
% Remove last 13 data points because we know they are messed up
RESraw(end-13:end,:) = [];
FMraw(end-13:end,:) = [];


% Read in the actual points tested and convert to format that can be simulated
testPoints_V = [RESraw(:,8), RESraw(:,7), RESraw(:,3), RESraw(:,2), RESraw(:,1)];   % actual tested points in units: [s(V), w(deg), P1(V), P2(V), P3(V)]
testPoints = [(testPoints_V(:,1)-4.077)*(1/0.033418), testPoints_V(:,2)-125, testPoints_V(:,3:5) * (testParams.TRmax/10)];   % convert it to units: [mm, deg, Pa, Pa, Pa]

% Calculate the predicted forces at each point in testPoints
FMpred = calcModelFM(testPoints, params);

% Separate points into separate configurations
cfg1 = zeros(size(testPoints(1,:)));    % seed with row of zeros
cfg2 = zeros(size(testPoints(1,:)));    % seed with row of zeros
cfg3 = zeros(size(testPoints(1,:)));    % seed with row of zeros
cfg4 = zeros(size(testPoints(1,:)));    % seed with row of zeros
FM1raw = zeros(1,2); FM2raw = zeros(1,2); FM3raw = zeros(1,2); FM4raw = zeros(1,2);
FM1pred = zeros(1,2); FM2pred = zeros(1,2); FM3pred = zeros(1,2); FM4pred = zeros(1,2);
for j = 1:length(testPoints)
    if abs(testPoints(j,1:2) - [0, 0]) < [2, 2]
        cfg1 = [cfg1; testPoints(j,:)];
        FM1raw = [FM1raw; FMraw(j,:)];
        FM1pred = [FM1pred; FMpred(j,:)]; 
    elseif abs(testPoints(j,1:2) - [5, 0]) < [2, 2]
        cfg2 = [cfg2; testPoints(j,:)];
        FM2raw = [FM2raw; FMraw(j,:)];
        FM2pred = [FM2pred; FMpred(j,:)]; 
    elseif abs(testPoints(j,1:2) - [0, 20]) < [2, 2]
        cfg3 = [cfg3; testPoints(j,:)];
        FM3raw = [FM3raw; FMraw(j,:)];
        FM3pred = [FM3pred; FMpred(j,:)]; 
    elseif abs(testPoints(j,1:2) - [5, 20]) < [2, 2]
        cfg4 = [cfg4; testPoints(j,:)];
        FM4raw = [FM4raw; FMraw(j,:)];
        FM4pred = [FM4pred; FMpred(j,:)]; 
    else
        disp('ERROR something wasnt classified!')
        disp(testPoints(j,:))
    end    
end
cfg1 = cfg1(2:end, :);    % remove row 1 seed of zeros
cfg2 = cfg2(2:end, :);    % remove row 1 seed of zeros
cfg3 = cfg3(2:end, :);    % remove row 1 seed of zeros
cfg4 = cfg4(2:end, :);    % remove row 1 seed of zeros
FM1raw = FM1raw(2:end, :); FM2raw = FM2raw(2:end, :); FM3raw = FM3raw(2:end, :); FM4raw = FM4raw(2:end, :);
FM1pred = FM1pred(2:end, :); FM2pred = FM2pred(2:end, :); FM3pred = FM3pred(2:end, :); FM4pred = FM4pred(2:end, :);


% Remove offset manually, using averages of measured offset
FM1 = -1 * ( FM1raw - FMraw(1,:) );
FM2 = -1 * ( FM2raw - [16.4, -0.022] - FMraw(1,:) );
FM3 = -1 * ( FM3raw - [-1.1, -0.021] - FMraw(1,:) );
FM4 = -1 * ( FM4raw - [16.2, -0.041] - FMraw(1,:) );
FM = [FM1; FM2; FM3; FM4];


% Calulate total error at each point
error = ( FM - FMpred );
sqerror = diag(error * error');
squerror1 = error(:,1).^2;
squerror2 = error(:,2).^2;


RMSE(1,1) = sqrt(sum(squerror1(1 : length(FM1(:,1)))) / (length(FM1(:,1))));
RMSE(1,2) = sqrt(sum(squerror2(1 : length(FM1(:,1)))) / (length(FM1(:,1))));
RMSE(2,1) = sqrt(sum(squerror1(length(FM1(:,1))+1 : length(FM1(:,1))+length(FM2(:,1)))) / (length(FM2(:,1))));
RMSE(2,2) = sqrt(sum(squerror2(length(FM1(:,1))+1 : length(FM1(:,1))+length(FM2(:,1)))) / (length(FM2(:,1))));
RMSE(3,1) = sqrt(sum(squerror1(length([FM1;FM2])+1 : length([FM1;FM2;FM3]))) / (length(FM3(:,1))));
RMSE(3,2) = sqrt(sum(squerror2(length([FM1;FM2])+1 : length([FM1;FM2;FM3]))) / (length(FM3(:,1))));
RMSE(4,1) = sqrt(sum(squerror1(length([FM1;FM2;FM3])+1 : end)) / (length(FM4(:,1))));
RMSE(4,2) = sqrt(sum(squerror2(length([FM1;FM2;FM3])+1 : end)) / (length(FM4(:,1))));
end