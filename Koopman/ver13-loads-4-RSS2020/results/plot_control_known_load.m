% plot_control_known_load
%   Creates plot of trajectory following for a known load from mpcData
%   files.
%   
%   Must first select mpcData files from file explorer (can select several 
%   simultaneously and it will display the different ones on a legend),
%   it generates a plot of that data, and you can save a cell array
%   containing all of (the paths to) the mpcData files you used.

% select data file(s)
[ datafile_name , datafile_path ] = uigetfile( '*.mat' , 'Choose data file(s) for merging...' , 'multiselect' , 'on' );

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

%% Generate plot

% set axis limits
xbounds = [-150,150];
ybounds = [-150,150];
zbounds = [-700,-400];

figure;

subplot(3,3,[1,2,4,5,7,8]);
hold on;
plot3( mpcData{1}.R(3:end,1) , mpcData{1}.R(3:end,3) , mpcData{1}.R(3:end,2) , 'LineWidth' , 3 );
plot3( mpcData{1}.Y(3:end,7) , mpcData{1}.Y(3:end,9) , mpcData{1}.Y(3:end,8) , 'LineWidth' , 3 );
hold off;
grid on; box on;
view(-37.5,15);
% axis equal;
xlim(xbounds);
ylim(ybounds);
zlim(zbounds);
set(gca,'DataAspectRatio',[1 1 1])
legend( { 'Reference' , 'Actual' } , 'Location' , 'southoutside' )

subplot(3,3,3); % xy-projection
hold on;
plot( mpcData{1}.R(3:end,1) , mpcData{1}.R(3:end,3) , 'LineWidth' , 3 );
plot( mpcData{1}.Y(3:end,7) , mpcData{1}.Y(3:end,9) , 'LineWidth' , 3 );
hold off;
grid on; box on;
xlim(xbounds);
ylim(ybounds);
% axis equal;
set(gca,'DataAspectRatio',[1 1 1])
title('xy-projection');

subplot(3,3,6); % xz-projection
hold on;
plot( mpcData{1}.R(3:end,1) , mpcData{1}.R(3:end,2) , 'LineWidth' , 3 );
plot( mpcData{1}.Y(3:end,7) , mpcData{1}.Y(3:end,8) , 'LineWidth' , 3 );
hold off;
grid on; box on;
xlim(xbounds);
ylim(zbounds);
% axis equal;
set(gca,'DataAspectRatio',[1 1 1])
title('xz-projection');

subplot(3,3,9); %yz-projection
hold on;
plot( mpcData{1}.R(3:end,3) , mpcData{1}.R(3:end,2) , 'LineWidth' , 3 );
plot( mpcData{1}.Y(3:end,9) , mpcData{1}.Y(3:end,8) , 'LineWidth' , 3 );
hold off;
grid on; box on;
xlim(ybounds);
ylim(zbounds);
% axis equal;
set(gca,'DataAspectRatio',[1 1 1])
title('yz-projection');
