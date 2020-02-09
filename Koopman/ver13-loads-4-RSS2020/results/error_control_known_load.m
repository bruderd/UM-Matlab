% error_control_known_load
%   Calculates RMSE error for trajectory following for a known load from mpcData
%   files.
%   
%   Must first select mpcData files from file explorer (can select several 
%   simultaneously and it will display the different ones on a legend),
%   it generates a plot of that data, and you can save a cell array
%   containing all of (the paths to) the mpcData files you used.

% select data file(s)
[ datafile_name , datafile_path ] = uigetfile( 'C:\Users\danie\Google Drive\PhD\Research\Labview\RSS2020\Matlab\*.mat' , 'Choose data file(s)...' , 'multiselect' , 'on' );

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

%% Caclulate RMSE error over each trial

% set the reference trajectory and start/stop indices
ref = mpcData{1}.R;
zorigin = ref(3,:);
Tref = 20; % time trajectories should take. EDIT FOR FASTER TRAJECTORIES!!
kmin = 3;   % index of first values that arent trash
kmax = floor( Tref / mpcData{1}.model.params.Ts );  % number of steps to get to the end of the reference trajectory (assuming it never finishes late)

RMSE = zeros( 2 , length( mpcData ) );
for i = 1 : length( mpcData )
    kmax_i = min( kmax , mpcData{i}.K(end) );
    RMSE(1,i) = mpcData{i}.Wreal(kmin,end);    % put load in first row
    terror = sqrt( sum( ( mpcData{i}.Y(kmin:kmax_i,7:9) - ref(kmin:kmax_i,:) ).^2 , 2 ) );
    RMSE(2,i) = sqrt( sum( terror.^2 ) / length(terror) ); % RMSE over whole trial
    RMSE(3,i) = mean( terror );     % average error over whole trial
end