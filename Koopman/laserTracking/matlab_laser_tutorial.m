% Matlab laser tracking tutorial

addpath('C:\Program Files\MATLAB\R2017a\toolbox\imaq');

%% Access and configure a device.
vid = videoinput('winvideo', 2, 'MJPG_320x240');   % '1' is built in webcam
vid.FramesPerTrigger = 1;
vid.TriggerRepeat = Inf;
triggerconfig(vid,'manual')

%% Create the laser figure window.
laserFig = figure;
hBox = plot([0 0 1 1 0], [0 1 1 0 0], 'b-');
hold on

% Set up calibration screen. Modify the cursor so it does not
% interfere with the calibration.
hTarget = plot(0, 0, 'yo');
ax = gca;
ax.Color = [0, 0, 0];
laserFig.Color = [0, 0, 0];
% laserFig.Menubar = 'none';
laserFig.DoubleBuffer = 'on';
laserFig.Pointer = 'custom';
laserFig.PointerShapeCData = NaN(16, 16);

%% Display positioning information.
posText = sprintf('%s\n%s', ...
    'Position the camera and ensure the blue box', ...
    'is the only thing in the camera''s view.');
infoText = text(0, -0.2,  posText, 'Color', [1 1 1]);
axis([-0.2 1.2 -0.2 1.2])
axis equal

% Using the preview window, request that the camera be positioned such
% that the view is of the blue box and little else.
preview(vid)
smallFigPos = laserFig.Position;
laserFig.Position = get(0, 'ScreenSize');
disp('Waiting for camera to be positioned...press any key to continue.')
pause

%% Perform Calibration
% Start the acquisition and create a new figure to display
% calibration results in a MATLAB SPY plot.
start(vid)
spyFig = figure;

% Target 1...
figure(laserFig);
hTarget.XData = 0;
hTarget.YData = 0;
sound(1), pause(2)
sound(1), trigger(vid);
acqResults{1} = getdata(vid, 1);

[xCalib(1), yCalib(1), laserSights] = util_findlaser(acqResults{1});
figure(spyFig);
spy(laserSights)
title('Target 1: Suspected Laser Sighting')

% Target 2...
figure(laserFig);
hTarget.XData = 0;
hTarget.YData = 1;
sound(1), pause(2)
sound(1), trigger(vid);
acqResults{2} = getdata(vid, 1);

[xCalib(2), yCalib(2), laserSights] = util_findlaser(acqResults{2});
figure(spyFig);
spy(laserSights)
title('Target 2: Suspected Laser Sighting')

% Target 3...
figure(laserFig);
hTarget.XData = 1;
hTarget.YData = 1;
sound(1), pause(2)
sound(1), trigger(vid);
acqResults{3} = getdata(vid, 1);

[xCalib(3), yCalib(3), laserSights] = util_findlaser(acqResults{3});
figure(spyFig);
spy(laserSights)
title('Target 3: Suspected Laser Sighting')

% Target 4...
figure(laserFig);
hTarget.XData = 1;
hTarget.YData = 0;
sound(1), pause(2)
sound(1), trigger(vid);
acqResults{4} = getdata(vid, 1);

[xCalib(4), yCalib(4), laserSights] = util_findlaser(acqResults{4});
figure(spyFig);
spy(laserSights)
title('Target 4: Suspected Laser Sighting')

% Close the SPY plot and stop the acquisition.
close(spyFig)
stop(vid);

%% Show calibration result
% Target 1 results...
calibFig = figure;
util_plotpos(acqResults{1}, xCalib(1), yCalib(1));

% Target 2 results...
util_plotpos(acqResults{2}, xCalib(2), yCalib(2));

% Target 3 results...
util_plotpos(acqResults{3}, xCalib(3), yCalib(3));

% Target 4 results...
util_plotpos(acqResults{4}, xCalib(4), yCalib(4));

% Close the figure illustrating calibration results.
close(calibFig)

%% Do laser tracking
% Update instructions on laser screen.
figure(laserFig);
infoText.String = 'Move the laser pointer within the blue box.';

% Start the acquisition. For each iteration:
%
% * output a sound to indicate a frame is about to be acquired
% * trigger the device
% * process the acquired image and locate the laser
% * convert pixel coordinates to MATLAB axis coordinates
laser.x = [];
laser.y = [];
start(vid)
for i = 1:100
    % Acquire an image frame and determine the
    % camera pixel coordinates.
    sound(1), trigger(vid);
    frame = getdata(vid, 1);
    [x, y] = util_findlaser(frame);

    if ~isnan(x) && ~isnan(y)
        % If coordinates were valid, ensure the camera pixel coordinate
        % was in the calibration range.
        x = max([x min(xCalib([1 2]))]);
        x = min([x max(xCalib([3 4]))]);
        y = min([y max(yCalib([1 4]))]);
        y = max([y min(yCalib([2 3]))]);

        % Determine spatial transformation from the unit square calibration points.
        tform = cp2tform([xCalib(:) yCalib(:)], [0 0; 0 1; 1 1; 1 0], 'projective');
        xyScreen = tformfwd([x, y], tform);
        xScreen = xyScreen(1);
        yScreen = xyScreen(2);

        % Ensure the new coordinates remain within the unit square.
        xScreen = min([xScreen 1]);
        xScreen = max([xScreen 0]);
        yScreen = min([yScreen 1]);
        yScreen = max([yScreen 0]);

        % Store the new MATLAB axis coordinates.
        laser.x = [laser.x(:); xScreen];
        laser.y = [laser.y(:); yScreen];
    end
end

% Plot the tracked laser positions.
laserFig.Position = smallFigPos;
plot(laser.x, laser.y, 'r*');

%% Close the laser figure.
input('Input anything to clear the laser figure and vid:')

close(laserFig);

% Stop the acquisition, remove the object from memory,
% and clear the variable.
stop(vid)
delete(vid)
clear vid

rmpath('C:\Program Files\MATLAB\R2017a\toolbox\imaq')


