% Load in test data and do stuff with it
clear
load('testData_1.mat')

FMraw = SyncedForceTorque.data;
FM = zeros(size(FMraw));
FM(:,1) = FMraw(:,1) + 2.02;
FM(:,2) = FMraw(:,2) - 0.101;