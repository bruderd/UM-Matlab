% main_laserTracking

%% Set Parameters

% set of calibration points?

%% Initialize Camera

% Access and configure a device.
vid = videoinput('winvideo', 2, 'MJPG_320x240');   % '1' is built in webcam
vid.FramesPerTrigger = 1;
vid.TriggerRepeat = Inf;
triggerconfig(vid,'manual')



%% Calibrate

[xCalib, yCalib] = calibrate(vid);      % perform cali

%% Run Tracking

%% Plot Results

%% Save Results (optional)

