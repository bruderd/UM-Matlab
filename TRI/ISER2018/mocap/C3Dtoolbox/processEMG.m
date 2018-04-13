% Process raw EMG values
% 
% EMG(t) -> e(t)
% 
% Tim Dorn
% June 2009
% 
% --------------------------------------------------------------------
% Usage: eVec = processEMG(EMGVec, {ProcessTask1, Value1, ProcessTask2, Value2, ...} )
% --------------------------------------------------------------------
% 
% Inputs:   EMGVec.data = Raw EMG data from extractRawEMG.m
%           EMGVec.time = time vector (sec)
%           EMGVec.name = string of muscle label
%           ProcessTaskX = Processing task X
%           ValueX = Value for processing task X
%           C3Dkey* = if this is given, the event lines will be added
%                     to the processed EMG plot
% 
% Outputs:  eVec = processed EMG structure (containing .time & .data & .name)
% 
% 
% Notes
% -----
% 
% Processing tasks are performed in the order they are given:
% 
% Task = 'REMDC' = Remove DC offsets        Value = []
% Task = 'RECT' = Full wave rectification   Value = []
% Task = 'REMDCRECT' = Remove DC offset & 
%                 full wave rectification   Value = []
% Task = 'HPF' = High pass filter           Value = [FilterOrder, Freq(Hz)]
% Task = 'LPF' = Low pass filter            Value = [FilterOrder, Freq(Hz)]
% Task = 'BPF' = Band pass filter           Value = [FilterOrder, FreqLow(Hz), FreqHigh(Hz)]
% Task = 'TKE' = TKE filter                 Value = []
% Task = 'NORM' = MVC Normalization         Value = MVC (uV)
% Task = 'NORM1' = Normalization to 1       Value = []
% Task = 'ABOVEZERO' = Force EMG > 0        Value = []
% Task = 'SQRT' = Squareroot EMG            Value = []
% Task = 'MULTIPLY' = Multiply EMG          Value = MultipleNumber
% Task = 'REMOVESPIKES' = Remove Spikes     Value = []
% Task = 'PLOT' = Plot2Screen value         Value = 0(OFF), 1(ON defualt), 
%                                                   'C3DKEY' variable string 
%                                                   (ON - SEPERATE EVENT INTERVALS 
%                                                         ON DIFFERENT SUBPLOTS)
% Task = 'SAVE' = Save plots to emf & fig   Value = [indices of process ops to not plot]
% Task = 'HIDE' = Hide line indices         Value = [](OFF default), (ON)
% Task = 'VERTLINES' = Plot event lines     Value = C3Dkey (default = [])
%               only active when PLOT = 1
% 
% 
% Note that the HPF and LPF options, a zero-phase Butterworth filter is used
% --------------------------------------------------------------------

function eVec = processEMG(EMGVec, processingTasks)

usage = 'Usage: eVec = processEMG(EMGVec, {ProcessTask1, Value1, ProcessTask2, Value2, ... } )';


                        
% Set up some initial parameters and do some initial checks
% ---------------------------------------------------------
larg = length(processingTasks);
    
if nargin < 2 || mod(larg, 2) == 1,
    disp(usage)
    return
end

plot2screen = 1;        % default values
saveplots = 0;
hide = [];
ploteventKey = [];



% Collate tasks
% -------------
numTasks = larg/2;
for i=1:numTasks
   tasks{i,1} = processingTasks{2*(i-1)+1};    % task name
   tasks{i,2} = processingTasks{2*(i-1)+2};    % task value
end
eVec = EMGVec;
legtxt = {'Raw EMG'};
glob = [];
key = [];
v=0;
% fprintf('\nMuscle: %s\n', upper(EMGVec.name));


% Perform tasks
% -------------
for i=1:numTasks
    
    switch upper(tasks{i,1})
        case 'REMDC'
            dispTXT = sprintf('Task %d: Remove DC Offset\n', i);
            eVec.data = removeDC(eVec);
            legtxt = [legtxt ; {'Removed DC'}];
            
            
        case 'RECT'
            dispTXT = sprintf('Task %d: Full Wave Rectification\n', i);
            eVec.data = fullWaveRect(eVec);
            legtxt = [legtxt ; {'Rectified'}];
            
            
         case 'REMDCRECT'
            dispTXT = sprintf('Task %d: Remove DC Offset & Full Wave Rectification\n', i);
            eVec.data = removeDCfullWaveRect(eVec); 
            legtxt = [legtxt ; {'RemDC & Rectified'}];
                
            
        case 'HPF'
            order = tasks{i,2}(1);
            freq = tasks{i,2}(2);
            dispTXT = sprintf('Task %d: High Pass Filter @ %d Hz (Order %d Zero Phase BW Filter)\n', ...
                i, freq, order);
            eVec.data = myfilter('high', eVec, order, freq);
            legtxt = [legtxt ; {sprintf('HPF @ %dHz (Ord=%d)', freq, order)}];
            
            
        case 'LPF'
            order = tasks{i,2}(1);
            freq = tasks{i,2}(2);
            dispTXT = sprintf('Task %d: Low Pass Filter @ %d Hz (Order %d Zero Phase BW Filter)\n', ...
                i, freq, order);
            eVec.data = myfilter('low', eVec, order, freq);
            legtxt = [legtxt ; {sprintf('LPF @ %dHz (Ord=%d)', freq, order)}];
            
            
        case 'BPF'
            order = tasks{i,2}(1);
            freqs = [tasks{i,2}(2), tasks{i,2}(3)];
            dispTXT = sprintf('Task %d: Band Pass Filter @ %d-%d Hz (Order %d Zero Phase BW Filter)\n', ...
                i, freqs(1), freqs(2), order);
            eVec.data = myfilter('bandpass', eVec, order, freqs);
            legtxt = [legtxt ; {sprintf('BPF @ %d-%d Hz (Ord=%d)', ...
                freqs(1), freqs(2), order)}];
            
            
        case 'TKE'
            dispTXT = sprintf('Task %d: TKE Filter\n', i);
            eVec.data = TKEfilter(eVec); 
            legtxt = [legtxt ; {'TKE filtered'}];
            
            
        case 'NORM'
            MVC = tasks{i,2}(1);
            dispTXT = sprintf('Task %d: MVC Normalization (Max MVC = %.3f uV)\n', ...
                i, MVC);
            eVec.data = MVCnorm(eVec, MVC);
            if max(eVec.data) > 1
                fprintf('WARNING: Normalized data [%s] > 1...\n', EMGVec.name);
            else
%                 fprintf('GOOD!: Normalized data [%s] < 1...\n', EMGVec.name);
            end
            
            if min(eVec.data) < 0
                fprintf('WARNING: Normalized data [%s] < 0...\n', EMGVec.name);
            else
%                 fprintf('GOOD!: Normalized data [%s] > 0...\n', EMGVec.name);
            end
            
            % Shape between 0 and 1
            eVec.data = min(1, eVec.data);
            eVec.data = max(0, eVec.data);
            
            legtxt = [legtxt ; {'Normalized'}];
            
            
        case 'NORM1'
            dispTXT = sprintf('Task %d: Normalization to 1\n', i);
            eVec.data = eVec.data / max(eVec.data);
            
            
        case 'ABOVEZERO'
            eVec.data = max(0, eVec.data);
            legtxt = [legtxt ; {'EMG Above Zero'}];
            
            
        case 'SQRT'
            eVec.data = sqrt(eVec.data);
            legtxt = [legtxt ; {'SQRT'}];    
            
            
        case 'MULTIPLY'
            f = tasks{i,2};
            eVec.data = f * eVec.data;
            legtxt = [legtxt ; {sprintf('MULTIPLY BY %.3f', f)}];       
            
            
        case 'REMOVESPIKES'
            eVec.data = removeEMGSpikes(eVec.data);
            legtxt = [legtxt ; {'Spikes Removed'}];
            
            
        case 'PLOT'
            plot2screen = tasks{i,2};
            v=v+1;
            
            
        case 'SAVE'
            saveplots = tasks{i,2};
            v=v+1;
            
            
        case 'HIDE'
            hide = tasks{i,2};
            v=v+1;    
            
            
        case 'VERTLINES'
            ploteventKey = tasks{i,2};
            v=v+1;
            
            
        otherwise
            error('Error: UNKNOWN TASK (%s)...\n', upper(tasks{i,1}));
    end
    glob(:,i-v) = eVec.data';
%     fprintf(dispTXT);

end



% Produce outputs
% ---------------
if isstruct(plot2screen)
    key = plot2screen;
    plot2screen = 2;
end

if plot2screen > 0 
    titl = sprintf('%s EMG: %s', eVec.c3dFile, eVec.name);
    figure('Name', titl, 'NumberTitle', 'off');
    hold on
    
    if plot2screen == 1
        h1 = plot(EMGVec.time, EMGVec.data, 'k', 'LineWidth', 1);
        h2 = plot(eVec.time, glob, 'LineWidth', 1);
        h = [h1;h2];
        if length(hide) > length(h)
            error('Trying to hide more objects than plotting (or not in range)\n')
        end
        delete(h(hide))
        xlim([EMGVec.time(1) EMGVec.time(end)])
        if ~isempty(ploteventKey)
            plotEventLines(ploteventKey, 'times0', 0);
        end
        legtxt = removeIndices(legtxt, hide);
        legend(legtxt, 'Location', 'Best')
        hline(0, 'k--');
        xlabel('time (sec)')
        title(titl, 'FontSize', 12)
        plotbrowser
    end
    
    if plot2screen == 2
        n = key.interval.numInterval;
        for i = 1:n
            ind = key.interval.Aframe0(i,1) : key.interval.Aframe0(i,2);
            t{i} = EMGVec.time(ind);
            rawCell{i} = EMGVec.data(ind);
            globCell{i} = glob(ind,:);
            hs = subplot(n+1,1,i);
            if i == key.interval4Time
                set(hs, 'Color', [0.99 0.99 0.97]);
            end
            hold on
            h1 = plot(t{i}, rawCell{i}, 'k', 'LineWidth', 1);
            h2 = plot(t{i}, globCell{i}, 'LineWidth', 1);
            h = [h1;h2];
            if length(hide) > length(h)
                error('Trying to hide more objects than plotting (or not in range)\n')
            end
            delete(h(hide))
            xlim([EMGVec.time(ind(1)) EMGVec.time(ind(end))])
            hline(0, 'k--');
            xlabel('time (sec)')
            ylabel(sprintf('Interval #%d', i))
            legtxt = removeIndices(legtxt, hide);
            legend(legtxt, 'Location', 'Best')
            legend hide
            if i == 1
                title(sprintf('%s             (%s -> %s)', titl, key.event.txt{1}, key.event.txt{1}), 'FontSize', 12)
            end
            plotbrowser
        end
        
        % calculate and plot average values
        tAvr = t{key.interval4Time};
        lt = length(tAvr);
        nOp = size(globCell{1},2);       % number of emg processing operations
        rawAvr = [];
        globAvr = [];
        for i = 1:n
            rawSamp(:, i) = resamp2(rawCell{i}', lt)';
            globSamp(:, nOp*(i-1)+1:nOp*i) = resamp2(globCell{i}', lt)';
        end
        rawAvr = mean(rawSamp,2);
        for i = 1:nOp
            globAvr(:,i) = mean(globSamp(:, i:nOp:n*nOp),2);
        end
        
        subplot(n+1,1,n+1)
        hold on
        h1 = plot(tAvr, rawAvr, 'k', 'LineWidth', 1);
        h2 = plot(tAvr, globAvr, 'LineWidth', 1);
        h = [h1;h2];
        if length(hide) > length(h)
            error('Trying to hide more objects than plotting (or not in range)\n')
        end
        delete(h(hide))
        xlim([tAvr(1) tAvr(end)])
        hline(0, 'k--');
        xlabel('time (sec)')
        legend(legtxt, 'Location', 'Best')
        legend hide
        ylabel('AVERAGE', 'FontSize', 12)
        
        % Overwrite eVec output with AVERAGE values
        eVec.SplitStrides.timeAvr1 = tAvr;              % time
        eVec.SplitStrides.timeAvr2 = tAvr - tAvr(1);    % time starting at 0
        eVec.SplitStrides.rawAvr = rawAvr;
        eVec.SplitStrides.processedAvr = globAvr(:,end);
    end
    
    if saveplots
        saveas(gcf, sprintf('%s_%s.emf', eVec.c3dFile, eVec.name), 'emf');
        saveas(gcf, sprintf('%s_%s', eVec.c3dFile, eVec.name), 'fig');
        print('-dpsc2', sprintf('%s_EMG.ps', eVec.c3dFile), '-append', sprintf('-f%s', num2str(gcf)))
    end
end








% ---------------------------------
% EMG PROCESSING TASK SUB-FUNCTIONS
% ---------------------------------

function data = removeDC(eVec)
    data = detrend(eVec.data);

    
function data = fullWaveRect(eVec)
    data = abs(eVec.data);   
    
    
function data = removeDCfullWaveRect(eVec)
    data = abs(detrend(eVec.data));   
   
    
function data = TKEfilter(eVec)
    x = eVec.data;
    l = length(x);
    
    for n = 1:l
        if n == 1
            data(n,1) = x(n)^2 - x(n+1)*x(n);
        elseif n == l
            data(n,1) = x(n)^2 - x(n)*x(n-1);
        else
            data(n,1) = x(n)^2 - x(n+1)*x(n-1);
        end
    end
 
    
function data = myfilter(filterType, eVec, order, freq)
    maxF = (1 / (eVec.time(2) - eVec.time(1))) / 2;
    [b_l,a_l] = butter(order, freq/maxF, filterType);
    data = filtfilt(b_l, a_l, eVec.data);  

    
function data = MVCnorm(eVec, MVC)
    data = eVec.data / MVC; 

