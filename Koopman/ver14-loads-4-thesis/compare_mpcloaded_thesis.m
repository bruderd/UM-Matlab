% compare_mpcloaded_thesis
%   Simulates a bilinear loaded and nonloaded model for the planar arm
%   doing a trajectory following task.
%
%   This is not a general purpose function.

%% Choose whether to save siulation results or not
saveon = false;

%% Load the models

% load the Arm system class
temp_arm = load('C:\Users\danie\OneDrive\Documents\MATLAB\UM-Matlab\Koopman\ver14-loads-4-thesis\systems\thesis-arm-markers_gravload1_3-mods_1-links_20hz\thesis-arm-markers_gravload1_3-mods_1-links_20hz.mat');
Arm = temp_arm.Arm;

% load the model classes
% temp_loaded = load('C:\Users\danie\OneDrive\Documents\MATLAB\UM-Matlab\Koopman\ver14-loads-4-thesis\systems\thesis-arm-markers_gravload1_3-mods_1-links_20hz\models\bilinear_poly-3_n-6_m-3_del-0_2020-06-21_16-39.mat');
temp_loaded = load('C:\Users\danie\OneDrive\Documents\MATLAB\UM-Matlab\Koopman\ver14-loads-4-thesis\systems\thesis-arm-markers_grav-endload-01_3-mods_1-links_20hz\models\bilinear_poly-3_n-6_m-3_del-0_2020-06-21_23-31.mat');
temp = load('C:\Users\danie\OneDrive\Documents\MATLAB\UM-Matlab\Koopman\ver14-loads-4-thesis\systems\thesis-arm-markers_noload_3-mods_1-links_20hz\models\bilinear_poly-3_n-6_m-3_del-0_2020-06-09_16-43.mat');
Ksysid_loaded = temp_loaded.sysid_class;  % bilinear model with load incorporated
Ksysid = temp.sysid_class;  % bilinear model without load incorporated

%% Construct an MPC controller for each model

% loaded bilinear mpc controller
Kmpc_loaded = Kmpc( Ksysid_loaded ,...
        'horizon' , 10 ,...
        'input_bounds' , 0.9*[-pi,pi],... % DO NOT USE INF
        'input_slopeConst' , [0.1e-1],... %1e-1 ,...
        'input_smoothConst' , [],... %[1e-1] ,...
        'state_bounds' , [] ,...
        'cost_running' , 10 ,...   % 0.1
        'cost_terminal' , 100 ,...  % 100
        'cost_input' , 0.1 * [ 3e-2 , 2e-2 , 1e-2 ]' ,...    % 1e-1
        'projmtx' , Ksysid_loaded.model.C(end-1:end,:) ,...  % just end effector
        'load_obs_horizon' , 20 ,...   % only needed for loaded models
        'load_obs_period' , 10 );   % only needed for loaded models
Ksim_loaded = Ksim( Arm , Kmpc_loaded );

% bilinear mpc controller (same as the one in chapter 4)
Kmpc = Kmpc( Ksysid ,...
        'horizon' , 10 ,...
        'input_bounds' , [],... % DO NOT USE INF
        'input_slopeConst' , [0.1e-1],... %1e-1 ,...
        'input_smoothConst' , [],... %[1e-1] ,...
        'state_bounds' , [] ,...
        'cost_running' , 10 ,...   % 0.1
        'cost_terminal' , 100 ,...  % 100
        'cost_input' , 0.1 * [ 3e-2 , 2e-2 , 1e-2 ]' ,...    % 1e-1
        'projmtx' , Ksysid.model.C(end-1:end,:) );% ,...  % just end effector
Ksim = Ksim( Arm , Kmpc );

%% Load in reference trajectory

% Circle
temp_ref = load('C:\Users\danie\OneDrive\Documents\MATLAB\UM-Matlab\Koopman\ver14-loads-4-thesis\trajectories\files\circle_c0-0p7_r0p3_15sec.mat');
ref = temp_ref.ref;


%% Sumulate each controller performing task over range of loads

loads = -pi/3 : pi/3 : pi/3;
res_loaded = cell( size(loads) );
res = cell( size(loads) );
for i = 1 : length( loads )
    res_loaded{i} =  Ksim_loaded.run_trial_mpc( ref.y , [0 0 0 0 0 0] , [0 0 0] , [ 1 , loads(i) ]);
    res{i} =  Ksim.run_trial_mpc( ref.y , [0 0 0 0 0 0] , [0 0 0] , [ 1 , loads(i) ]);
end

%% save results (OPTIONAL)
if saveon
    sim_folder = [ 'systems' , filesep , Kmpc_lin.model.params.sysParams.sysName , filesep , 'simulations' ];
    mkdir( sim_folder , ref.name );
    filepath = [ sim_folder , filesep , ref.name ];
    
    filename_loaded = [ filepath , filesep , Kmpc_loaded.params.classname , '.mat' ];
    save( filename_loaded , 'res_loaded' );
    
    filename = [ filepath , filesep , Kmpc.params.classname , '.mat' ];
    save( filename , 'res' );
end


%% Plot results

% colormap
colormap lines;
cmap = colormap;

% axis limits
axis_limits = [ -50 , 50 , 20 , 120 ];
x_tick_marks = [ -50 , -25 , 0 , 25 , 50 ];
y_tick_marks = [ 0 , 25 , 50 , 75 , 100 ];

% Grid points for gravity direction arrows
arrow_len = 10;
[ x_grid , y_grid ] = meshgrid( axis_limits(1) : arrow_len : axis_limits(2) , axis_limits(3) : arrow_len : axis_limits(4) );

% Just the trajectories
figure; % plot in cm instead of meters
for i = 1 : length( loads )
    % direction of gravity
    u_grid = -ones( size(x_grid) ) * arrow_len * sin( res{i}.W(end,2) );
    v_grid = ones( size(y_grid) ) * arrow_len * cos( res{i}.W(end,2) );
    
    subplot(1,length(loads),i);
    hold on;
    quiver( x_grid , y_grid , u_grid , v_grid , 'Color' , [0.65 0.65 0.65] );
    plot( 100 * res_loaded{i}.R(:,1) , 100 * res_loaded{i}.R(:,2) , 'Color' , [0 0 0 0.5] , 'LineWidth' , 2 );    % reference trajectory
    plot( 100 * res_loaded{i}.Y(:,end-1) , 100 * res_loaded{i}.Y(:,end) , ':' , 'Color' , [ cmap(4,:) , 1 ] , 'LineWidth' , 2);    % system trajectory
    plot( 100 * res{i}.Y(:,end-1) , 100 * res{i}.Y(:,end) , 'Color' , [ cmap(2,:) , 1 ] , 'LineWidth' , 2 );    % non-loaded system trajectory
    hold off;
    axis( axis_limits );
    xticks( x_tick_marks );
    yticks( y_tick_marks );
    set(gca, 'YDir','reverse')
    daspect([1 1 1]);   % make axis ratio 1:1
    box on; grid on;
    xlabel('$\hat{\alpha}$ (cm)' , 'Interpreter' , 'Latex');
    if i == 1
        ylabel('$\hat{\beta}$ (cm)' , 'Interpreter' , 'Latex');
    end
    title(['$w = [' , num2str(res{i}.W(end,1)) , ',' , num2str(res{i}.W(end,2)) , ']^\top $' ] , 'Interpreter' , 'Latex');
end

%% Plot trajectories and error

figure; % plot in cm instead of meters
for i = 1 : length( loads )
    % direction of gravity
    u_grid = -ones( size(x_grid) ) * arrow_len * sin( res{i}.W(end,2) );
    v_grid = ones( size(y_grid) ) * arrow_len * cos( res{i}.W(end,2) );
    
    % trajectories
    subplot(2,length(loads),i);
    hold on;
    quiver( x_grid , y_grid , u_grid , v_grid , 'Color' , [0.65 0.65 0.65] );
    plot( 100 * res_loaded{i}.R(:,1) , 100 * res_loaded{i}.R(:,2) , 'Color' , [0 0 0 0.5] , 'LineWidth' , 2 );    % reference trajectory
    plot( 100 * res_loaded{i}.Y(:,end-1) , 100 * res_loaded{i}.Y(:,end) , 'Color' , [ cmap(4,:) , 1 ] , 'LineWidth' , 2);    % system trajectory
    plot( 100 * res{i}.Y(:,end-1) , 100 * res{i}.Y(:,end) , 'Color' , [ cmap(2,:) , 1 ] , 'LineWidth' , 2 );    % non-loaded system trajectory
    hold off;
    axis( axis_limits );
    xticks( x_tick_marks );
    yticks( y_tick_marks );
    set(gca, 'YDir','reverse')
    daspect([1 1 1]);   % make axis ratio 1:1
    box on; grid on;
    xlabel('$\hat{\alpha}$ (cm)' , 'Interpreter' , 'Latex');
    if i == 1
        ylabel('$\hat{\beta}$ (cm)' , 'Interpreter' , 'Latex');
    end
    title(['$w = [' , num2str(res{i}.W(end,1)) , ',' , num2str(res{i}.W(end,2)) , ']^\top $' ] , 'Interpreter' , 'Latex');
    
    % error
    subplot(2,length(loads),i+length(loads));
    hold on;
    plot( res_loaded{i}.T(1:end-1) , 100 * res_loaded{i}.err , 'Color' , [ cmap(4,:) , 1 ] , 'LineWidth' , 2);    % system trajectory
    plot( res{i}.T(1:end-1) , 100 * res{i}.err , 'Color' , [ cmap(2,:) , 1 ] , 'LineWidth' , 2);    % non-loaded system trajectory
    hold off;
    box on; grid on;
    daspect([1 5 1]);   % make axis ratio 1:1
    ylim([0,50]);
    yticks([0,25,50]);
    xlabel('Time (seconds)');
    if i == 1
        ylabel('Tracking Error (cm)');
    end
end






