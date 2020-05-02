function err = evaluate_nonlinear_models( data )
%evaluate_nonlinear_models: Evaluates the performance of koopman model for soft
%robot system verses Matlab sysid toolbox models and neural network

% If data argument not supplied, prompt uset to load data from file
if ~exist('data', 'var')
    % Prompt user to identify data file
    [data_file,data_path] = uigetfile('..\dataFiles\*.mat');
    all_data = load([data_path, data_file]); % must be a .mat file
    data = all_data.val;    % isolate just the validation data trials
end

%% simulate nonlinear Koopman model

% load in the nonlinear koopman model
matcontents = load('n-2_m-3_del-0_2020-04-29_21-20.mat');
ksysid = matcontents.sysid_class;
Ts = ksysid.params.Ts;   % timestep
koop = cell( size(data) );
for i = 1 : length(data)
    koop{i}.t = data{i}.t;
    koop{i}.u = ksysid.scaledown.u( data{i}.u );
    ykoop0 = ksysid.scaledown.y( data{i}.y(1,:) );
    [ ~ , koop{i}.y ] = ode45( @(t,y) ksysid.model.F_func( y , koop{i}.u( floor(t/Ts)+1 ,:)' ) , koop{i}.t , ykoop0' );
    
    % scale data back up
    koop{i}.u = data{i}.u;
    koop{i}.y = ksysid.scaleup.y( koop{i}.y );
    
    % distance error at each timestep
    koop{i}.err = sqrt( sum( ( data{i}.y - koop{i}.y ).^2 , 2 ) );
end


%% simulate neural network model

nnet = cell( size(data) );
for i = 1 : length(data)
    nnet{i}.t = data{i}.t;
    nnet{i}.u = data{i}.u;
    [ ynnet , ~ ] = nnet_model_nio( nnet{i}.u' , [0 0 0]' );
    nnet{i}.y = ynnet';
    
    % distance error at each timestep
    nnet{i}.err = sqrt( sum( ( data{i}.y - nnet{i}.y ).^2 , 2 ) );
end


%% simulate NLARX model

% load in the NLARX model from file
load('nlarx1.mat');
nlarx = cell( size(data) );
for i = 1 : length( data )
    nlarx{i}.t = data{i}.t;
    nlarx{i}.u = data{i}.u;
    udata = iddata( [] , nlarx{i}.u , nlarx1.Ts ); 
    ydata = sim( nlarx1 , udata );
    nlarx{i}.y = ydata.y;
    
    % distance error at each timestep
    nlarx{i}.err = sqrt( sum( ( data{i}.y - nlarx{i}.y ).^2 , 2 ) );
end

%% simulate Hamm-Weiner model

% load in the NLARX model from file
load('nlhw1.mat');
nlhw = cell( size(data) );
for i = 1 : length( data )
    nlhw{i}.t = data{i}.t;
    nlhw{i}.u = data{i}.u;
    udata = iddata( [] , nlhw{i}.u , nlhw1.Ts ); 
    ydata = sim( nlhw1 , udata );
    nlhw{i}.y = ydata.y;
    
    % distance error at each timestep
    nlhw{i}.err = sqrt( sum( ( data{i}.y - nlhw{i}.y ).^2 , 2 ) );
end

%% evaluate the total average error for each model

% sum the error from all the trials
koop_allerr = [];
nnet_allerr = [];
nlarx_allerr = [];
nlhw_allerr = [];
% koop_toterr = 0;
% nnet_toterr = 0;
% nlarx_toterr = 0;
% nlhw_toterr = 0;
% totsteps = 0;
for i = 1 : length( data )
    % concatenate all the error vectors
    koop_allerr = [ koop_allerr ; koop{i}.err ];
    nnet_allerr = [ nnet_allerr ; nnet{i}.err ];
    nlarx_allerr = [ nlarx_allerr ; nlarx{i}.err ];
    nlhw_allerr = [ nlhw_allerr ; nlhw{i}.err ];
    
%     koop_toterr = koop_toterr + sum( koop{i}.err );
%     nnet_toterr = nnet_toterr + sum( nnet{i}.err );
%     nlarx_toterr = nlarx_toterr + sum( nlarx{i}.err );
%     nlhw_toterr = nlhw_toterr + sum( nlhw{i}.err );
    
%     totsteps = totsteps + length( data{i}.t );
end

% calculate average error for each model
totsteps = size( koop_allerr , 1 );
err.koop.ave = sum( koop_allerr ) / totsteps;
err.nnet.ave = sum( nnet_allerr ) / totsteps;
err.nlarx.ave = sum( nlarx_allerr ) / totsteps;
err.nlhw.ave = sum( nlhw_allerr ) / totsteps;
% koop_averr = koop_toterr / totsteps;
% nnet_averr = nnet_toterr / totsteps;
% nlarx_averr = nlarx_toterr / totsteps;
% nlhw_averr = nlhw_toterr / totsteps;

% calculate standard deviation
err.koop.std = std( koop_allerr );
err.nnet.std = std( nnet_allerr );
err.nlarx.std = std( nlarx_allerr );
err.nlhw.std = std( nlhw_allerr );



%% plot the trajectory comparison for each validation trial


%% plot the error comparison over all validation trials

% Define custom colors:
cb(1,:) = [69,117,180]/255;
cb(2,:) = [145,191,219]/255;
cb(3,:) = [254,224,144]/255;
cb(4,:) = [252,141,89]/255;
cb(5,:) = [215,48,39]/255;


% bar chart showing TOTAL NRMSE across all states and all trials
figure
hold on
bar(1 , err.koop.ave * 2.54 , 'FaceColor' , cb(1,:) );
bar(2 , err.nnet.ave * 2.54 , 'FaceColor' , cb(2,:) );
bar(3 , err.nlarx.ave * 2.54 , 'FaceColor' , cb(3,:) );
bar(4 , err.nlhw.ave * 2.54 , 'FaceColor' , cb(4,:) );
xticks(1:4);
xtickangle(45);
xticklabels( {'Koopman', 'Neural Network', 'NLARX', 'Ham.-Wiener'} );
ylabel('Average Error (cm)');
% ylim([0 10]);
% eb = errorbar(1:numSystems-1, [mean_NRMSE_alltrials(1:3)' * 100, mean_NRMSE_alltrials(5:end)' * 100], [std_NRMSE_alltrials(1:3)' * 100, std_NRMSE_alltrials(5:end)' * 100], '.', 'CapSize', 14, 'LineWidth', 1.5, 'Color', 'k');
eb = errorbar(1:4, [err.koop.ave , err.nnet.ave , err.nlarx.ave , err.nlhw.ave] * 2.54, [err.koop.std , err.nnet.std , err.nlarx.std , err.nlhw.std] * 2.54, '.', 'CapSize', 14, 'LineWidth', 1.5, 'Color', 'k');
box on
hold off
