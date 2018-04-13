% Extract, filter, and plot data from a saved OpenSim data file (*.mot or *.sto)
% into a Matlab structure
% Tim Dorn
% July 2009
% 
% --------------------------------------------------------------------
% Usage: [out, plots2make, labels, dataFile, hp, xl] = extractMotFile(ProcessTask1, Value1, ProcessTask2, Value2, ...)
% --------------------------------------------------------------------
% 
% Inputs:   ProcessTaskX = Processing task X
% 
%       Task: 'FILE'        Value = filename string (default = choose file manually)
%       Task: 'XAXIS'       Value = field for x axis (default = 'time')
%       Task: 'FILT'        Value = low pass filter frequency (default = 0)
%       Task: 'PLOT'        Value = plot variables
%                                -1 = off
%                                 0 = user selectable (default)
%                                 [n1 n2] = plot the columns specified by the array
%       Task: 'FIG'         Value = figure to plot into (default = 123)
%       Task: 'VERTLINES'   Value = C3D event key (from getEvents.m)
%       Task: 'EMG'         Value = EMG mot file (from batchEMGprocess.m)
%       Task: 'GRAPHMIN'    Value = Minimum value for vertical axis of graphs
%       Task: 'GRAPHMAX'    Value = Maximum value for vertical axis of graphs
%       Task: 'HLINE'       Value = Vector of horizontal lines for each index
%       Task: 'FIBRELEN'    Value = OpenSim model structure (from Load_OSIM)
%       Task: 'ALLON1FIG'   Value = 0 = seperate figures, 1 = one figure (default)
% 
% 
% Outputs:  out.name:     trial name
%           out.labels:   extracted data labels
%           out.data:     extracted data (after filtering)
%           out.filtFreq: low pass filter frequency
%           performPlots: indices of plots extracted
%           labels:       all labels from dataFile
%           dataFile:     filename used
%           hp:           plot handle
%           
% 
% Notes:    If no input arguments are given i.e. data = extractMotFile, 
%           the function allows you to select a file for plotting and will
%           superimpose over existing plots if a common variable is being
%           plotted.
% 
% --------------------------------------------------------------------

function [out, plots2make, labels, dataFile, hp, xl] = extractMotFile(varargin)

larg = length(varargin);
usage = 'Usage: [out, plots2make, labels, dataFile] = extractMotFile(ProcessTask1, Value1, ProcessTask2, Value2, ...)';

if mod(larg, 2) == 1,
    disp(usage)
    return
end


% Collate tasks
% -------------
numTasks = length(varargin)/2;
for i=1:numTasks
   tasks{i,1} = varargin{2*(i-1)+1};    % task name
   tasks{i,2} = varargin{2*(i-1)+2};    % task value
end


% Set defaults
% ------------
file = [];
dataFile = [];
dataFilePreLoaded = [];
xaxis = 'time';
yaxis = '';
performPlots = 0;
freq = 0;
figNum = 123;
C3Dkey = [];
emgFile = [];
ALLON1FIG = 1;
hl = [];
model_addfibrelen = [];
penn = 0;
supT = [];
maxPlot = 0;
figPrep = 1;
dataInput = [];
saveFile = 0;
addexcitation = [];
addactivation = [];
addkinematics = [];
model_addtendonstrain = [];
writeTitle = 1;
LineW = '';
LineC = '';
LineS = '';
hline0 = 1;
model_maxfom = [];
calcmean = 0;
mult = 1;
add = 0;
differentiate = 0;
addcontrolconstraints = [];
normModel = [];
addwork = 0;
timecrop = [];



% Set Variables
% -------------
for i=1:numTasks
    
    switch upper(tasks{i,1})
        case 'FILE'
            dataFile = tasks{i,2};

        case 'FILEPRELOADED'
            dataFilePreLoaded = tasks{i,2};            
            
        case 'XAXIS'
            xaxis = tasks{i,2};
            
        case 'YAXIS'
            yaxis = tasks{i,2};            
            
        case 'FILT'
            freq = tasks{i,2};

        case 'PLOT'
            performPlots = tasks{i,2};

        case 'FIG'
            figNum = tasks{i,2};
            
        case 'VERTLINES'
            C3Dkey = tasks{i,2};
            
        case 'EMG'
            emgFile = tasks{i,2};
            EMGdataset = extractMotFile('file', emgFile, 'plot', -1);
            timeEMG = getData('time', EMGdataset);
            
        case 'GRAPHMIN'
            graphMin = tasks{i,2};

        case 'GRAPHMAX'
            graphMax = tasks{i,2};
            
        case 'HLINE'
            hl = tasks{i,2};
            
        case 'FIBRELEN'
            model_addfibrelen = tasks{i,2};
            
        case 'PENNATION'
            penn = tasks{i,2};     
            
        case 'ALLON1FIG'
            ALLON1FIG = tasks{i,2};

        case 'SUPERTITLE'
            supT = tasks{i,2};
            
        case 'MAX'
            maxPlot = tasks{i,2};    
            
        case 'FIGPREP'
            figPrep = tasks{i,2};    
            
        case 'DATAINPUT'
            dataInput = tasks{i,2};      
            
        case 'SAVEFILE'
            saveFile = tasks{i,2}; 
            
        case 'ADDEXCITATION'
            addexcitation = tasks{i,2}; 
            excLab = extractMotFile('file', addexcitation, 'plot', -1);
            
        case 'ADDACTIVATION'
            addactivation = tasks{i,2}; 
            actLab = extractMotFile('file', addactivation, 'plot', -1);

        case 'ADDTENDONSTRAIN'
            model_addtendonstrain = tasks{i,2};
            
        case 'ADDCONTROLCONSTRAINTS'
            addcontrolconstraints = tasks{i,2};
            if ~isempty(addcontrolconstraints)
                cc = xml_read(addcontrolconstraints);
            end
            
        case 'ADDKINEMATICS'
            addkinematics = tasks{i,2}{1};
            Qind = tasks{i,2}{2};
            
            
        case 'WRITETITLE'
            writeTitle = tasks{i,2}; 

        case 'LINEW'
            LineW = tasks{i,2}; 

        case 'LINEC'
            LineC = tasks{i,2}; 

        case 'LINES'
            LineS = tasks{i,2}; 
            
        case 'HLINE0'
            hline0 = tasks{i,2}; 
            
        case 'MAXFOM'
            model_maxfom = tasks{i,2};
            
        case 'CALCMEAN'
            calcmean = tasks{i,2};
 
        case 'MULTIPLY'
            mult = tasks{i,2};
            
        case 'ADD'
            add = tasks{i,2};    
            
        case 'DIFF'
            differentiate = tasks{i,2};
            
        case 'NORMALIZE'
            normModel = tasks{i,2}{1};
            normField = tasks{i,2}{2};

        case 'ADDWORK'
            addwork = tasks{i,2};    
            
        case 'TIMECROP'
            timecrop = tasks{i,2};
            
        otherwise
            error('Error: UNKNOWN TASK (%s)...\n', upper(tasks{i,1}));
    end
end

if freq < 0, 
    error('Error: filtering frequency must be > 0');
end

if isempty(dataInput) && isempty(dataFilePreLoaded)
    
    if isempty(dataFile)
       [file, path] = uigetfile({'*.mot;*.sto', 'OpenSim Motion Files (*.mot, *.sto)'; ...
           '*.*',  'All Files (*.*)'}, 'Select Motion File');
        dataFile = [path, file];
    else
        file = dataFile;
    end



    % Extract motion data from external file
    % --------------------------------------
    nlrows = 1;
    nhead = 1;

    try
        % rename nans to be detectable. REQUIRES ssr.exe in system32 directory
%         eval(sprintf('!ssr 0 "1.#QNAN000000000000000" "NaN" "%s"', dataFile));
%         eval(sprintf('!ssr 0 "1.#IND0000000000000000" "NaN" "%s"', dataFile));
    catch
        fprintf('NaNs in the file will cause an error...\n')
    end

    fid = fopen(dataFile, 'r');
    if fid < 0
        fprintf('\nERROR: %s could not be opened for reading...\n\n', dataFile);
        return
    end

    buffer = fgetl(fid);
    while strcmpi(strtrim(buffer), 'endheader') == 0,
        if nhead == 1,
            out.name = buffer;
        end
        nhead = nhead + 1;
        if ~isempty(findstr(buffer, 'nColumns')),
            % get number of columns from the header
            ncols = strread(buffer, '%*s%d', 'delimiter', '=');
        end
        if ~isempty(findstr(buffer, 'nRows')),
            % get number of rows from the header
            nrows = strread(buffer, '%*s%d', 'delimiter', '=');
        end
        buffer = fgetl(fid);
    end
    fclose(fid);

    [labels, x, y] = readColData(dataFile, ncols, nhead, nlrows);
    out.labels = labels;
    out.data = [x, y];
else
    out = dataInput;
end

if ~isempty(dataFilePreLoaded)
   out = dataFilePreLoaded;
   labels = out.labels;
end


% perform multiplication/addition operations here
out.data = [out.data(:,1) out.data(:,2:end) * mult + add];
xData = getData(xaxis, out);
sampFreq = floor(1 / (xData(2) - xData(1)));



% Filter Data (optional)
% ----------------------

if freq > 0,
    fprintf('Filtering Data @ %d Hz\n', freq);
    for i = 1:ncols
        if ~strcmpi(out.labels{i}, 'time'),
            out.data(:,i) = smooth(out.data(:,i), freq, sampFreq);
        end
    end
end
out.freq = freq;



% Plot extracted data
% -------------------
if ~isempty(performPlots) && performPlots(1) >=0, 
    
%     ALLON1FIG = 1;          % 0 = seperate figures, 1 = one figure
    
    % Range of plot colours
    col = [0 0 1;
           1 0 0;
           0 0 0;
           0.2 0.5 0.1;
           0.5 0.4 0.68;
           1 0.7 0.2;
           0.7 0.7 0.3
           0.31 0.50 0.67
           0.23 0.62 0.62
           0.75 0.50 0.41
           0.39 0.31 0.10];
       
    lw = [2 2 3];                             % Range of line widths
    ls = [{'-'}, {'--'}, {':'}];     % ['.', 'x', '*', '+', 'o',];      % Range of line styles
    
    for i = 1:length(out.labels)
        listLabels{i} = sprintf('[%d]  %s', i, out.labels{i});
        
        if isstruct(normModel)
            try
                if strcmpi(normField, 'Vmax')
                    out.data(:,i) = out.data(:,i) / ...
                        (normModel.Muscles.(strtok(out.labels{i}, '.')).(normField) * ...
                        normModel.Muscles.(strtok(out.labels{i}, '.')).optimal_fiber_length);
                else
                    out.data(:,i) = out.data(:,i) / ...
                        normModel.Muscles.(strtok(out.labels{i}, '.')).(normField);
                end
            catch
            end
        end

    end
    
    if performPlots(1) == 0,
        [plots2make,v] = listdlg('PromptString', 'Variables to plot:',...
                      'SelectionMode', 'multiple',...
                      'ListString', listLabels);
        if v == 0, return; end      % cancel button pressed
    else
        plots2make = performPlots;
    end
    
    l = length(plots2make);
    pl = ceil(sqrt(l));
                    
    % Main plotting loop
    for i = 1:l
        j = plots2make(i);

        if differentiate
            pp = spline(xData, out.data(:,j));
            ppd = diffpp(pp);
            out.data(:,j) = ppval(ppd, xData);
        end
        
        if figPrep == 1
            if ALLON1FIG == 0
                pl = 1;
                figure(100+j)
                
            elseif ALLON1FIG == 1
                h = figure(figNum);
                set(h,'PaperType','a4', 'PaperPositionMode','manual',...
                    'PaperOrientation','landscape', 'PaperUnits','centimeters',...
                    'PaperPosition', [0,0,29,20]);
                hs = subplot(pl, pl, i);
                hold on
                
            elseif length(ALLON1FIG) == 2
                h = figure(figNum);
                hs = subplot(ALLON1FIG(1), ALLON1FIG(1), ALLON1FIG(1)*(ALLON1FIG(2)-1)+i);
                hold on
            end
        end
        
        if maxPlot ~= 0
            maximizecustom(maxPlot);
        end
        
        hold on
%         , 'Name', sprintf('%s', out.labels{j}), 'NumberTitle', 'off');
        [a,b,c,d] = legend;
        
        if ~isempty(LineW)
            lw = LineW;
        end
        if ~isempty(LineC)
            col = LineC;
        end
        if ~isempty(LineS)
            ls = LineS;
        end
        
        hp = plot(xData, out.data(:,j), ...
            'Color', col(1+mod(length(d), size(col,1)), :), ...
            'LineWidth', lw(1+mod(length(d), length(lw))), ...
            'LineStyle', ls{1+mod(length(d), length(ls))});
        uistack(hp, 'top');
            
%     col{1+mod(length(d), length(col))})
        
        avg = mean(out.data(:,j));
        avgabs = mean(abs(out.data(:,j)));
        maxval = max(out.data(:,j));
        meanrms = rms(out.data(:,j));
        
%             leg = sprintf('%s (avg = %.2f, avgabs = %.2f, max = %.2f)', ...
%                 file, avg, avgabs, maxval);

        if calcmean
            hline(meanrms, 'r--', sprintf('RMS = %.3f\n', meanrms));
        end

        % The legend value denotes the legend in the matlab plotbrowser
        leg = out.labels{j};
        legend([d';{leg}], 'Location', 'Best', 'Interpreter', 'none');
        legend hide
        xlabel(xaxis, 'FontSize', 8, 'Interpreter', 'none')
        
        % plot yaxis only on left hand subplots
        if mod(i-1, pl) == 0
            ylabel(yaxis, 'FontSize', 8, 'Interpreter', 'none')
        end
        
        axis tight
        xl = xlim;
        
        if hline0
            hline(0, 'k:');
        end
        
        if ~isempty(hl)
            if length(hl) == 1
                hline(hl, 'k:');
            else
                hline(hl(j), 'k:');
            end
        end
        

        if isstruct(model_maxfom)
            fom = model_maxfom.Muscles.(strtok(labels{j}, '.')).max_isometric_force;
            hline(fom, 'k:');
        end
        
        
        axis tight
        croptime(timecrop, xl)
        yl = ylim;
        diff = 0.1*(yl(2) - yl(1));
        newyl = [yl(1)-diff, yl(2)+diff];
        ylim(newyl);

        if exist('graphMin', 'var')
            ylim([graphMin, yl(2)]);
        end
        
        if exist('graphMax', 'var')
            ylim([yl(1), graphMax]);
        end        
        

        
        
        % add kinematics
        if ~isempty(addkinematics)
            hold on
            xData2 = addkinematics.data(:, getIndex(addkinematics.labels, 'time'));
            yData2 = addkinematics.data(:, Qind);
            Ucol = [.3 .8 .3];

            [haxes, hline1, hline2] = plotyy(xData, out.data(:,j), xData2, yData2, ...
                'plot', 'plot');
            set(hline1, 'LineStyle', 'none')
            set(hline1, 'LineWidth', 2)
            set(hline2, 'LineWidth', 2, 'Color', Ucol)

            hold on
            out.haxes = haxes;

            axis tight
            croptime(timecrop, xl)
            ylabel('MomArm / JntAng', 'FontSize', 8, 'Interpreter', 'none')
            set(gca, 'YTickLabelMode', 'auto')
            set(gca, 'YTickMode', 'auto')
            ylim(newyl);
            
            axes(haxes(2));
            axis tight
            croptime(timecrop, xl)
            set(gca, 'YTickLabelMode', 'auto')
            set(gca, 'YTickMode', 'auto')
            yl2 = ylim;
            diff2 = 0.1*(yl2(2) - yl2(1));
            newyl2 = [yl2(1)-diff2, yl2(2)+diff2];
            ylim(newyl2);
        end
        
        
        
        
        
        % add excitation / activation values
        if ~isempty(addexcitation)
            hold on
            xData2 = excLab.data(:, getIndex(excLab.labels, 'time'));
            try
                yData2 = excLab.data(:, getIndex(excLab.labels, [strtok(out.labels{j}, '.') '.excitation']));
                Ucol = [.3 .8 .3];
            catch
                try
                    yData2 = excLab.data(:, getIndex(excLab.labels, strtok(out.labels{j}, '.')));
                    Ucol = [1 .64 .57];
                catch
                    yData2 = excLab.data(:, getIndex(excLab.labels, [strtok(out.labels{j}, '.') '.activation']));
                    Ucol = [.4 .6 .7];
                end
            end
            
            [haxes, hline1, hline2] = plotyy(xData, out.data(:,j), xData2, yData2, ...
                'plot', 'plot');
            set(hline1, 'LineStyle', 'none')
            set(hline1, 'LineWidth', 2)
            set(hline2, 'LineWidth', 1, 'Color', Ucol)

            hold on
            out.haxes = haxes;

%                 axes(haxes(1));  
            axis tight
            croptime(timecrop, xl)
            if mod(i-1, pl) == 0
                ylabel('Force (N) / Excitation', 'FontSize', 8, 'Interpreter', 'none')
            end
            set(gca, 'YTickLabelMode', 'auto')
            set(gca, 'YTickMode', 'auto')
            yl = ylim;
            ylim([0, newyl(2)]);
%             ylim(newyl);
            
            axes(haxes(2));
            axis tight
            croptime(timecrop, xl)
            ylim([0, 1.1]);
            set(gca, 'YTickLabelMode', 'auto')
            set(gca, 'YTickMode', 'auto')
        end
        
        
        if ~isempty (addcontrolconstraints)
            % get muscle control constraints
            for k = 1:length(cc.objects.ControlLinear)
                t = []; minP = []; maxP = [];
                if strcmpi(cc.objects.ControlLinear(k).ATTRIBUTE.name, ...
                        sprintf('%s.excitation', strtok(out.labels{j}, '.')))
                    try 
                        for m = 1:length(cc.objects.ControlLinear(k).min_nodes.ControlLinearNode)
                            t(m) = cc.objects.ControlLinear(k).min_nodes.ControlLinearNode(m).t;
                            minP(m) = cc.objects.ControlLinear(k).min_nodes.ControlLinearNode(m).value;
                            maxP(m) = cc.objects.ControlLinear(k).max_nodes.ControlLinearNode(m).value;
                        end
                    catch

                    end
                    break
                end
            end

            t = [xl(1) t xl(2)];
            if isempty(minP)
                minP = [0.01 minP 0.01];
            else
                minP = [minP(1) minP minP(end)];
            end
            if isempty(maxP)
                maxP = [1 maxP 1];
            else
                maxP = [maxP(1) maxP maxP(end)];
            end               

            % plot shaded control constraints
            hold on
            [ha hb hc] = shadedplot(t, minP, maxP, [1 1 0.8]);
            set(hb, 'LineStyle', '--', 'LineWidth', 1, 'Color', [1 0.93 0.8])
            set(hc, 'LineStyle', '--', 'LineWidth', 1, 'Color', [1 0.93 0.8])
            uistack(ha(2), 'bottom')
            alpha(0.2)
        end
            
                    
        if ~isempty(addactivation)
            xData3 = actLab.data(:, getIndex(actLab.labels, 'time'));
            try
                yData3 = actLab.data(:, getIndex(actLab.labels, [strtok(out.labels{j}, '.') '.activation']));
                Acol = [.71 .87 .5];
            catch
                yData3 = actLab.data(:, getIndex(actLab.labels, strtok(out.labels{j}, '.')));
                Acol = [1 .64 .57];
            end
            hold on
            plot(xData3, yData3, 'LineStyle', '-', 'LineWidth', 1, 'Color', Acol)
        end       
        
        
        % add tendon strain values
        % (assumes the plotting of tendon length with model as second variable input)
        if ~isempty(model_addtendonstrain)
            hold on
            lst = model_addtendonstrain.Muscles.(strtok(out.labels{j}, '.')).tendon_slack_length;
            hline(lst, 'k:');
            tendonStrain = (out.data(:,j) - lst)./lst;
            if any(tendonStrain<0)
                fprintf('Warning: Tendon strain < 0 at some stage for %s\n', out.labels{j})  
            end
            tendonStrain(tendonStrain<0) = 0;
            
            [haxes, hline1, hline2] = plotyy(xData, out.data(:,j), xData, tendonStrain, ...
                'plot', 'plot');
            set(hline1, 'LineStyle', 'none')
            set(hline1, 'LineWidth', 2)
            set(hline2, 'LineWidth', 1)
%             set(hline1, 'LineStyle', 'none')       % remove tendon strain line
            
            hold on
            out.haxes = haxes;

%                 axes(haxes(1));  
            axis tight
            croptime(timecrop, xl)
            if mod(i-1, pl) == 0
                ylabel('Lengths (m) / TendonStrain', 'FontSize', 8, 'Interpreter', 'none')
            end
            set(gca, 'YTickLabelMode', 'auto')
            set(gca, 'YTickMode', 'auto')
            yl = ylim;
            ylim([0, 1.1*yl(2)]);
            
            % assume here we are also plotting muscle fiblen (plotVarArms(9) && plotVarArms(20))
            ofl = model_addtendonstrain.Muscles.(strtok(out.labels{j}, '.')).optimal_fiber_length;
            ShadePlotForEmpahsisHoriz([0.5*ofl, 1.5*ofl], 'y', 0.2);
            hline(ofl, 'k:');
            
            axes(haxes(2));
            axis tight
            croptime(timecrop, xl)
            ylim([0, 0.1]);
            set(gca, 'YTickLabelMode', 'auto')
            set(gca, 'YTickMode', 'auto')
        end

        
        % add work = integral of power
        % (assumes the plotting of power - this will overlay the integral)
        if addwork
            hold on
            power = cumtrapz(xData, out.data(:,j));  % trapezoidal integration of power

            [haxes, hline1, hline2] = plotyy(xData, out.data(:,j), xData, power, ...
                'plot', 'plot');
            set(hline1, 'LineStyle', 'none')
            set(hline1, 'LineWidth', 2)
            set(hline2, 'LineWidth', 1)
            
            hold on
            out.haxes = haxes;
            axis tight
            croptime(timecrop, xl)
            if mod(i-1, pl) == 0
                ylabel('Fibre Power/Work', 'FontSize', 8, 'Interpreter', 'none')
            end
            set(gca, 'YTickLabelMode', 'auto')
            set(gca, 'YTickMode', 'auto')
            ylim(newyl);

            axes(haxes(2));
            croptime(timecrop, xl)
            axis tight
            set(gca, 'YTickLabelMode', 'auto')
            set(gca, 'YTickMode', 'auto')
            yl2 = ylim;
            diff2 = 0.1*(yl2(2) - yl2(1));
            newyl2 = [yl2(1)-diff2, yl2(2)+diff2];
            ylim(newyl2);
        end        
        
        
        % plot fibre length / pennation angles
        % (assumes the plotting of muscle fibre length with model as second variable input)
        if ~isempty(model_addfibrelen)
            ofl = model_addfibrelen.Muscles.(strtok(out.labels{j}, '.')).optimal_fiber_length;
            ShadePlotForEmpahsisHoriz([0.5*ofl, 1.5*ofl], 'y', 0.2);
            hline(ofl, 'k:', sprintf('%f\n', ofl));
            graphMin = 0;
            ylim([graphMin, 1.5*ofl]);
            
            if penn
                hold on
                muscWidth = model_addfibrelen.Muscles.(out.labels{j}).optimal_fiber_length * ...
                    sin(model_addfibrelen.Muscles.(out.labels{j}).pennation_angle);
                pennAng = asin(muscWidth./out.data(:,j)) * 180/pi;
                
                if abs(sum(imag(pennAng))) > 0
                   fprintf('Warning: 90deg pennation angle for %s\n', out.labels{j})  
                end
                
                [haxes, hline1, hline2] = plotyy(xData, out.data(:,j), xData, ...
                    real(pennAng), 'plot', 'plot');
                set(hline1, 'LineStyle', 'none')
                set(hline2, 'LineStyle', '--')
                set(hline1, 'LineWidth', 2)
                set(hline2, 'LineWidth', 1)
                
                hold on
                ofl = model_addfibrelen.Muscles.(out.labels{j}).optimal_fiber_length;
%                 ShadePlotForEmpahsisHoriz([0.5*ofl, 1.5*ofl], 'y', 0.2);
%                 hline(ofl, 'b--', sprintf('%f\n', ofl));
                out.haxes = haxes;

%                 axes(haxes(1));  
                axis tight
                xlim(xl)
                if mod(i-1, pl) == 0
                    ylabel('FibLen (m) / PenAng (deg)', 'FontSize', 8, 'Interpreter', 'none')
                end
                set(gca, 'YTickLabelMode', 'auto')
                set(gca, 'YTickMode', 'auto')
                yl = ylim;
                ylim([0, 1.1*yl(2)]);
                axes(haxes(2));
                axis tight
                xlim(xl)                
                set(gca, 'YTickLabelMode', 'auto')
                set(gca, 'YTickMode', 'auto')
            end
        end


        

%             plotbrowser
        if writeTitle
            if ~isempty(C3Dkey)
                title(gca, sprintf('%s\n', strtok(out.labels{j}, '.')), ...
                    'Interpreter', 'none', 'FontSize', 10, 'FontWeight', 'Bold')
            else
                title(gca, strtok(out.labels{j}, '.'), ...
                    'Interpreter', 'none', 'FontSize', 10, 'FontWeight', 'Bold')
            end
        end

        % Check for EMG to superimpose...
        if ~isempty(emgFile)
            EMG_heightRatio = 0.3;
            EMGsingleMusc = getData(sprintf('%s.Processed', strtok(out.labels{j}, '.')), EMGdataset);
            if ~isempty(EMGsingleMusc)
%                 fprintf('Found %s in EMG dataset!\n', out.labels{j})
                maxEMG = max(EMGsingleMusc);
                GraphLimit = max(ylim) - min(ylim);
                maxEMGline = (GraphLimit*EMG_heightRatio) + min(ylim);
                scaledEMG = EMGsingleMusc .* (GraphLimit*EMG_heightRatio/maxEMG) + min(ylim);
                hold on
                hline(maxEMGline, 'c:');
                plot(timeEMG, scaledEMG, 'Color', [0.678 0.922 1.0])
            end
        end

        if saveFile && ~ALLON1FIG 
            saveas(gcf, sprintf('Figure%s.emf', num2str(gcf)), 'emf');
            saveas(gcf, sprintf('Figure%s.fig', num2str(gcf)), 'fig');
%           print('-dpsc2', sprintf('Figure%s.ps', num2str(gcf)), sprintf('-f%s', num2str(gcf)))
            close(100+j)
        end
    
        croptime(timecrop, xl)
        
        
        if ~isempty(C3Dkey)
            plotEventLines(C3Dkey, 'times0', 0);
        end
        
    end
end

if ~isempty(supT)
    suptitle(supT); 
end

if saveFile && ALLON1FIG 
    saveas(gcf, sprintf('Figure%s.emf', num2str(gcf)), 'emf');
    saveas(gcf, sprintf('Figure%s.fig', num2str(gcf)), 'fig');
%     print('-dpsc2', sprintf('Figure%s.ps', num2str(gcf)), sprintf('-f%s', num2str(gcf)))
end

% turn off legends for subplots (all except 1)
% if ALLON1FIG == 1
%     for i = 2:l
%         subplot(pl, pl, i)
%         legend off
%     end
% end




function croptime(timecrop, xl)

if ~isempty(timecrop)
    for k = 1:2
        if timecrop(k) < 0
            timecrop(k) = xl(k);
        end
    end
    xlim(timecrop);
else
    xlim(xl);
end


