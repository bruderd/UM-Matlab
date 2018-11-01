% main_laserTracking

%% Set Parameters

params = struct;

params.calPoints = [0 0; 0 1; 1 1; 1 0; 0.5 0.5];  % set of calibration points?
params.numFrames = 200;     % how long to perform tracking

%% Initialize Camera

% Access and configure a device.
vid = videoinput('winvideo', 2, 'MJPG_800x600');   % '1' is built in webcam
vid.FramesPerTrigger = 1;
vid.TriggerRepeat = Inf;
triggerconfig(vid,'manual')


%% Calibrate

calib = calibrate_laser(vid,params);      % perform cali

%% Run Tracking

laser = track_laser(vid, calib , params); % laser is struct with .x , .y , .t


% Stop the acquisition, remove the object from memory,
% and clear the variable.
delete(vid)
clear vid

%% Plot Results

resultPlot = figure;
hold on
plot(laser.x , laser.y, '*');    % draw connected plot of all the points the laser was recorded
plot(laser.x , laser. y);

%% Save Results (optional)

