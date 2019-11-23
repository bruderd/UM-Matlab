% example_KMPC.m

%%
%-------------------------------------------------------------------------%
%----------------- Identify Koopman Model from Data ----------------------%
%-------------------------------------------------------------------------%

% load in data from file
[ datafile_name , datafile_path ] = uigetfile( 'datafiles/*.mat' , 'Choose data file for sysid...' );
data4sysid = load( [datafile_path , datafile_name] );

% construct sysid class
sysid = sysid( data4sysid, ...
        'model_type' , 'linear' ,...
        'obs_type' , { 'poly' } ,...
        'obs_degree' , [ 3 ] ,...
        'snapshots' , Inf ,...
        'lasso' , [ 10 ] ,...
        'delays' , 0 ,...
        'isupdate' , false );

    
% train model(s)
sysid = sysid.train_models;

% validate model(s)
results = cell( size(sysid.candidates) );    % store results in a cell array
err = cell( size(sysid.candidates) );    % store error in a cell array 

if iscell(sysid.candidates)
    for i = 1 : length(sysid.candidates)
        [ results{i} , err{i} ] = sysid.valNplot_model( i );
    end
else
    [ results{1} , err{1} ] = sysid.valNplot_model;
end


%%
%-------------------------------------------------------------------------%
%----------------- Construct Koopman MPC Controller ----------------------%
%-------------------------------------------------------------------------%

mpc = mpc( sysid ,...
        'horizon' , 25 ,...
        'input_bounds' , [ 0 , 10 ],... 
        'input_slopeConst' , [0.5e-2],... 
        'input_smoothConst' , [1e-1],... 
        'state_bounds' , [] ,...
        'cost_running' , 0.1 ,...
        'cost_terminal' , 100 ,...
        'cost_input' , 0 ); 
    
%%
%-------------------------------------------------------------------------%
%-------------------------- Run Simulation -------------------------------%
%-------------------------------------------------------------------------%

% load in reference trajectory from file
[ reffile_name , reffile_path ] = uigetfile( 'ref-trajectories/*.mat' , 'Choose reference trajectory file...' );
temp = load( [reffile_path , reffile_name] );
ref = temp.ref;

% run simulation
y0 = [ 1 ,1 ];      % initial laser dot position
u0 = [ 0 , 0 , 0];  % initial input
results = mpc.run_simulation( ref , y0 , u0 );


%%
%-------------------------------------------------------------------------%
%---------------------------- Plot Results -------------------------------%
%-------------------------------------------------------------------------%

% Define color(s)
cb_blue = [55,126,184] ./ 255;  % blue

figure;
hold on;
plot( results.R(:,1) , results.R(:,2) , 'color' , [0.25 0.25 0.25 1] , 'LineWidth' , 2);
plot( results.Y(:,1) , results.Y(:,2) , 'color' , cb_blue , 'LineWidth' , 2);
xlabel('$x_1$ (cm)' , 'Interpreter' , 'Latex' , 'FontSize', 20);
ylabel('$x_2$ (cm)' , 'Interpreter' , 'Latex' , 'FontSize', 20);
hold off;
grid on; box on;





