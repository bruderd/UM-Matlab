% Matlab code to control myDaq valves
function controlValves(Pcontrol, sampleRate, timeBetweenPoints)

% Convert Pcontrol to TR control signal in Volts 

% Create a control signal to run at the designated sample rate



%% Run test
s = daq.createSession('ni');        % create a session

addAnalogInputChannel(s, 'myDAQ1', 0, 'Voltage');   % sets function of pin 0 of myDAQ1
addAnalogOutputChannel(s, 'myDaQ1', 1, 'Voltage');  % sets funtion of pin 1 of myDAQ1

s.Rate = 100;    % set scan rate to 100 samples per second
% s.DurationInSeconds = 10;   % run the aquisition for 10 seconds (only needed if we are only measuring)

queueOutputData(s, controlPressures);   % generate control signal

[captured_data, time] = s.startForeground();     % starts the aquisition and returns the results