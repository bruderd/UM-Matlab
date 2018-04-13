% Get the MVC reference voltages from a set of contractions
% Tim Dorn
% June 2008
% 
% --------------------------------------------------------------------
% Usage: [MVC, MVCset] = getMVC(C3Dkey, emgSetName, windowSize, *MVCmethod)
% --------------------------------------------------------------------
% 
% Inputs:   C3Dkey: the C3D key structure from getEvents
%           emgSetName: the label of EMG names contained in the EMG set
%               (this must be defined in loadlabels.m as glab.[emgSetName])
%           windowSize = gliding 'look ahead' window size (sec)
% 
%           *MVCmethod = calculating method of MVC (used once the maximum
%                       mean EMG time window has been found)
%               'MEAN' --> MVC = mean value (default value if not specified)
%               'RMS'  --> MVC = root mean square value
%               'MAX'  --> MVC = maximum value
% 
% 
% Outputs:  MVC = structure of MVC voltages (uV)    
%                   (MVC.[muscLabel] = value)
% 
%           MVCset = structure of cells containing information about the
%                    EMG labels
% 
% 
% Notes
% -----
% 
% EMG LABELS ARE CONTAINED IN: loadLabels.m
% 
% References: ABC or EMG (page 30)
% 
% Bolgla, L. A. and T. L. Uhl (2007). "Reliability of electromyographic
% normalization methods for evaluating the hip musculature." J Electromyogr
% Kinesiol 17(1): 102-11.
% 
% ---------------------------------------------------------------------

function MVC = getMVC(C3Dkey, emgSetName, windowSize, MVCmethod)

usage = 'Usage: MVC = getMVC(C3Dkey, emgSetName, windowSize, *MVCmethod)';


% Set up some initial parameters and do some initial checks
% ---------------------------------------------------------
if nargin == 3,
    MVCmethod = 'MEAN';
elseif nargin ~= 4,
    disp(usage)
    return
end



% Extract & Process MVC from EMG data
% -----------------------------------
[emg, MVCset] = extractRawEMG(C3Dkey, emgSetName);
numMuscles = length(MVCset);
MVC.info = sprintf('windowSize = %.2f sec, MVCmethod = %s', ...
    windowSize, MVCmethod);
MVC.array = [];
time = emg.(strtok(MVCset{1}{1}, '.()[]{}')).time;


figure('Name', 'MVC EMG Calculations', 'NumberTitle', 'off');
hold on 
for i = 1:numMuscles
    musc_name_orig = MVCset{i}{1};
    musc_name = strtok(musc_name_orig, '.()[]{}');
    musc = emg.(musc_name);
    subplot(numMuscles, 1, i)
    
    musc2 = processEMG(musc, ...
        'remdc', [], ...            % Remove DC Offset
        'hpf', [4, 20], ...         % High Pass Filter
        'rect', [], ...             % Rectification
        'plot', 0);                 % Don't Plot
    
    
    plot(time, musc2.data)
    ylabel('uV')
    xlim([time(1), time(end)])
    
    % Get the MVC value
    [MVCval, t1, t2] = getMVCvalue(musc2, windowSize, MVCmethod);
    
    hline(double(MVCval), 'r--');
    ShadePlotForEmpahsisVert([t1, t2], 'y', 0.3);
    MVC.(musc_name) = MVCval;
    MVC.array = [MVC.array, MVCval];
    title(sprintf('%s --> %s (MVC = %.4f)', ...
        musc_name_orig, MVCset{i}{2}, MVCval))
end
% fprintf('MVC Method Used: %s\n', upper(MVCmethod));






% ========================================================================
% SUBFUNCTION: getMVCvalue
% ========================================================================
% Steps through the EMG data with a defined windowSize (sec) and finds the
% window size(t1<t<t2) where the mean EMG amplitude is a maximum. 
% The MVC is then calculated using the MVCmethod.

function [mvc, t1, t2] = getMVCvalue(muscle, windowSize, MVCmethod)

timePeriod = muscle.time(2) - muscle.time(1);
frameWS = round(windowSize / timePeriod);       % frame window size
if frameWS > length(muscle.time),
    error('Error: Window size is greater than the entire trial...');
end

% Find time window corresponding to maximum mean EMG amplitude
% ------------------------------------------------------------
l = length(muscle.time) - frameWS;
mvc = 0;
for i = 1:l
    data = muscle.data(i:i+frameWS);
    mvctmp = mean(data);
    if mvctmp > mvc,
        mvc = mvctmp;
        dataWINDOW = data;
        t1 = muscle.time(i);
        t2 = muscle.time(i+frameWS);
    end
end



% Use MVC method to now calculate the MVC value
% ---------------------------------------------
switch upper(MVCmethod)
    
    case 'MEAN'
        mvc = mean(dataWINDOW);
        
    case 'RMS'
        mvc = rms(dataWINDOW);
        
    case 'MAX'
        mvc = max(dataWINDOW);
        
    otherwise
        error('Error: UNKNOWN MVCmethod (%s)...\n', upper(MVCmethod));
end

