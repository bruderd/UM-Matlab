% plot_control_unknown_load
%   Creates plot of trajectory following for an unknown load from mpcData
%   files.
%   
%   Must first select mpcData files from file explorer (can select several 
%   simultaneously and it will display the different ones on a legend),
%   it generates a plot of that data, and you can save a cell array
%   containing all of (the paths to) the mpcData files you used.

% select data file(s)
[ datafile_name , datafile_path ] = uigetfile( 'C:\Users\danie\Google Drive\PhD\Research\Labview\RSS2020\Matlab\mpc\dataFiles\*.mat' , 'Choose data file(s)...' , 'multiselect' , 'on' );

% load in the data files
min_trial_len = Inf;    % the length of the shortest trial loaded in (in steps not time)
if iscell( datafile_name )  % check if it's cell array
    data_merged = cell( 1 , length(datafile_name) );
    for i = 1 : length(datafile_name)
        temp = load( [datafile_path , datafile_name{i}] );
        mpcData{i} = temp.mpcData;
        min_trial_len = min( length( mpcData{i}.T ) , min_trial_len );
    end
else    % if not a cell array, turn it into 1x1 cell array
    data_merged = cell(1,1);
    temp = load( [datafile_path , datafile_name] );
    mpcData{1} = temp.mpcData;
    min_trial_len = length( mpcData{1}.T );
    disp('FYI, you only selected one file so your output cell array will have dimension 1.');
end

%% Generate plot

% set the reference trajectory (assumes all trials have same reference)
ref = mpcData{1}.R(1:min_trial_len,:);
zorigin = ref(3,:);

% initialze fields for legend
trials = gobjects( 3 , length(mpcData) + 1);
trialnames = cell( 1 , length(mpcData) + 1);
trialnames{1} = 'Reference';

% set axis limits
xbounds = [-150,150];
ybounds = [-150,150];
zbounds = [-150,150];

figure;
% set(gcf, 'Position',  [100, 100, 1400, 600])
colormap lines;
clr = colormap;
refclr = [1 1 1]*0.6;   % color for reference trajectory

% Load estimate verses time
subplot(2,2,[1,2])
hold on;
for i = 1 : length( mpcData )
    time = mpcData{i}.T(3:min_trial_len) - mpcData{i}.T(3);
    plot( time , mpcData{i}.Wreal(3,end) * ones(size(time)) , '--' , 'LineWidth' , 1.5 , 'Color' , clr(i,:) );
    plot( time , mpcData{i}.W(3:min_trial_len) , 'LineWidth' , 2 , 'Color' , clr(i,:));
%     trialnames{i+1} = ['Actual ($w =$ ', num2str(mpcData{i}.Wreal(3,end)) , ' g)'];    % change to W
end
grid on; box on;
hold off;
xlim([0,40]);
ylabel('Load Estimate, $\hat{w}$ (g)')
xlabel('Time (seconds)')

% Tracking error verses time
subplot(2,2,[3,4])
hold on;
for i = 1 : length( mpcData )
    time = mpcData{i}.T(3:min_trial_len) - mpcData{i}.T(3);
    terror = sqrt( sum( ( mpcData{i}.Y(3:min_trial_len,7:9) - ref(3:min_trial_len,:) ).^2 , 2 ) );
    plot( time , terror , 'LineWidth' , 2 , 'Color' , clr(i,:));
end
grid on; box on;
hold off;
ylim([0,100]);
xlim([0,40]);
ylabel('Tracking error (mm)')
xlabel('Time (seconds)')

% % 3d plot plus legend
% subplot(2,4,[3,4,7,8])
% hold on;
% trials(:,1) = plot3( ref(3:min_trial_len,1) , ref(3:min_trial_len,3) , ref(3:min_trial_len,2) - zorigin , 'LineWidth' , 3 , 'Color' , refclr );
% for i = 1 : length( mpcData )
%     trials(:,i+1) = plot3( mpcData{i}.Y(3:min_trial_len,7) , mpcData{i}.Y(3:min_trial_len,9) , mpcData{i}.Y(3:min_trial_len,8) - zorigin , 'LineWidth' , 2 , 'Color' , clr(i,:));
%     trialnames{i+1} = ['Actual ($w =$ ', num2str(mpcData{i}.Wreal(3,3)) , ' g)'];    % change to W
% end
% hold off;
% grid on; box on;
% view(-37.5+180,15);
% xlim(xbounds);
% ylim(ybounds);
% zlim(zbounds);
% set(gca,'DataAspectRatio',[1 1 1])
% legend( trials(1,:) , trialnames , 'Location' , 'northoutside' , 'Interpreter' , 'Latex' );








