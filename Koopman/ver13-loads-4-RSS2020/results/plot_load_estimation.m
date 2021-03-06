% plot_load_estimation
%   Creates a plot of the load estimates for the robot converging over time
%
%   Must first select mpcData files from file explorer (can select several 
%   simultaneously and it will display the different ones on a legend),
%   it generates a plot of that data, and you can save a cell array
%   containing all of (the paths to) the mpcData files you used.


% select data file(s)
[ datafile_name , datafile_path ] = uigetfile( 'C:\Users\danie\Google Drive\PhD\Research\Labview\RSS2020\Matlab\load_estimation\dataFiles\*.mat' , 'Choose data file(s)...' , 'multiselect' , 'on' );

% load in the data files
if iscell( datafile_name )  % check if it's cell array
    data_merged = cell( 1 , length(datafile_name) );
    for i = 1 : length(datafile_name)
        temp = load( [datafile_path , datafile_name{i}] );
        mpcData{i} = temp.mpcData;
    end
else    % if not a cell array, turn it into 1x1 cell array
    data_merged = cell(1,1);
    temp = load( [datafile_path , datafile_name] );
    mpcData{1} = temp.mpcData;
    disp('FYI, you only selected one file so your output cell array will have dimension 1.');
end

%% Reorder data in order of smallest to largest load

load = zeros( length(mpcData) , 2 );    % assumes 1-D load
temp = mpcData; % temporary copy of mpcData
for i = 1 : length( mpcData )
    load(i,:) = [ i , mpcData{i}.Wreal(3,end) ];
end
load_ordered = sortrows( load , 2 , 'ascend' );
for i = 1 : length( mpcData )
    mpcData{i} = temp{ load_ordered(i,1) };
end

%% Generate plot

figure;
% set(gcf, 'Position',  [100, 100, 1400, 600])
% set(0,'defaulttextinterpreter','latex')
colormap lines;
clr = colormap;
refclr = [1 1 1]*0.6;   % color for reference trajectory

% set time for trial
Tmax = 30; % time trajectories should take. EDIT FOR FASTER TRAJECTORIES!!
kmin = 3;   % index of first values that arent trash
kmax = floor( Tmax / mpcData{1}.model.params.Ts );  % number of steps to get to the end of the reference trajectory (assuming it never finishes late)

% initialze fields for legend
trials = gobjects( 1 , length(mpcData) );
trialnames = cell( 1 , length(mpcData) );

% set axis limits
xbounds = [0,Tmax];
ybounds = [0,300];

% Load estimate verses time
hold on;
for i = 1 : length( mpcData )
    time = mpcData{i}.T(kmin:kmax) - mpcData{i}.T(kmin);
    plot( time , mpcData{i}.Wreal(kmin,end) * ones(size(time)) , '--' , 'LineWidth' , 1.5 , 'Color' , clr(i,:) );
    trials(i) = plot( time , mpcData{i}.W(kmin:kmax) , 'LineWidth' , 2 , 'Color' , clr(i,:));
    trialnames{i} = ['\textsf{Actual load} $w =$ ', num2str(mpcData{i}.Wreal(kmin,end)) , ' g'];
end
grid on; box on;
hold off;
ylabel('Load Estimate, $\hat{w}$ (g)');
xlabel('Time (seconds)');
% legend( trials(1,:) , trialnames , 'Location' , 'northeast' , 'Interpreter' , 'Latex' );











