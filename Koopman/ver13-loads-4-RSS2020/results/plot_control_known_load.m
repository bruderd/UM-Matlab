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

% set the reference trajectory
ref = mpcData{1}.R;
zorigin = ref(3,:);
trialnames = cell( 1 , length(mpcData) + 1);
trialnames{1} = 'Reference';

% set axis limits
xbounds = [-150,150];
ybounds = [-150,150];
zbounds = [-150,150];

figure;

subplot(3,3,[1,2,4,5,7,8]);
hold on;
plot3( ref(3:end,1) , ref(3:end,3) , ref(3:end,2) - zorigin , 'LineWidth' , 3 );
for i = 1 : length( mpcData )
    plot3( mpcData{i}.Y(3:end,7) , mpcData{i}.Y(3:end,9) , mpcData{i}.Y(3:end,8) - zorigin , 'LineWidth' , 2 );
    trialnames{i+1} = ['Actual ($w =$ ', num2str(mpcData{i}.What(end,3)) , ' g)'];    % change to W
end
hold off;
grid on; box on;
view(-37.5,15);
xlim(xbounds);
ylim(ybounds);
zlim(zbounds);
set(gca,'DataAspectRatio',[1 1 1])
legend( trialnames , 'Location' , 'southoutside' , 'Interpreter' , 'Latex' );

subplot(3,3,3); % xy-projection
hold on;
plot( ref(3:end,1) , ref(3:end,3) , 'LineWidth' , 3 );
for i = 1 : length( mpcData )
    plot( mpcData{i}.Y(3:end,7) , mpcData{i}.Y(3:end,9) , 'LineWidth' , 2 );
end
hold off;
grid on; box on;
xlim(xbounds);
ylim(ybounds);
set(gca,'DataAspectRatio',[1 1 1])
title('xy-projection');

subplot(3,3,6); % xz-projection
hold on;
plot( ref(3:end,1) , ref(3:end,2) - zorigin , 'LineWidth' , 3 );
for i = 1 : length( mpcData )
    plot( mpcData{i}.Y(3:end,7) , mpcData{i}.Y(3:end,8) - zorigin , 'LineWidth' , 2 );
end
hold off;
grid on; box on;
xlim(xbounds);
ylim(zbounds);
set(gca,'DataAspectRatio',[1 1 1])
title('xz-projection');

subplot(3,3,9); %yz-projection
hold on;
plot( ref(3:end,3) , ref(3:end,2) - zorigin , 'LineWidth' , 3 );
for i = 1 : length( mpcData )
    plot( mpcData{i}.Y(3:end,9) , mpcData{i}.Y(3:end,8) - zorigin , 'LineWidth' , 2 );
end
hold off;
grid on; box on;
xlim(ybounds);
ylim(zbounds);
set(gca,'DataAspectRatio',[1 1 1])
title('yz-projection');
