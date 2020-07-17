function [ err , res ] = evaluate_all_models( data )
%evaluate_all_models: Evaluates the performance of koopman models for soft
%robot system verses Matlab sysid toolbox models and neural network

% If data argument not supplied, prompt uset to load data from file
if ~exist('data', 'var')
%     % Prompt user to identify data file
%     [data_file,data_path] = uigetfile('..\dataFiles\*.mat');
%     all_data = load([data_path, data_file]); % must be a .mat file
%     data = all_data.val;    % isolate just the validation data trials

    % default data
    temp = load('mean_noisedata.mat');
    data = temp.ndata;
end


%% simulate linear Koopman model

% load in the linear koopman model
matcontents = load('n-2_m-3_del-1_2020-04-30_19-13.mat');
ksysid_lin = matcontents.sysid_class;
Ts = ksysid_lin.params.Ts;   % timestep
kooplin = cell( size(data) );

% % (OPTIONAL) Resample the data coming in (should remove if not needed)
% for i = 1 : length( data )
%     tq = ( data{i}.t(1) : Ts : data{i}.t(end) )';
%     data{i}.u = interp1( data{i}.t , data{i}.u , tq );
%     data{i}.y = interp1( data{i}.t , data{i}.y , tq );
%     data{i}.t = tq;
% end

for i = 1 : length(data)
    kooplin{i}.t = data{i}.t;
    
    % scale down the data
    datasc.t = data{i}.t;
    datasc.u = ksysid_lin.scaledown.u( data{i}.u );
    datasc.y = ksysid_lin.scaledown.y( data{i}.y );
    
    [~,kooplin{i}.zeta] = ksysid_lin.get_zeta( datasc );
    kooplin{i}.z = zeros( length(kooplin{i}.t) , length( ksysid_lin.basis.full ) );   % preallocate  
    kooplin{i}.z(1,:) = ksysid_lin.lift.full( kooplin{i}.zeta(1,:)' )'; % initial value
    kooplin{i}.y = zeros( size( data{i}.y ) );  % preallocate
    kooplin{i}.y(1,:) = datasc.y(1,:);  % initial value
    
    % simulate the (scaled) system
    for j = 1 : length( kooplin{i}.t ) - 1
        kooplin{i}.z(j+1,:) = ( ksysid_lin.model.A * kooplin{i}.z(j,:)' + ksysid_lin.model.B * datasc.u(j,:)' )';
        kooplin{i}.y(j+1,:) = ( ksysid_lin.model.C * kooplin{i}.z(j+1,:)' )';
    end
    
    % scale data back up
    kooplin{i}.u = data{i}.u;
    kooplin{i}.y = ksysid_lin.scaleup.y( kooplin{i}.y );
    
    % distance error at each timestep
    kooplin{i}.err = sqrt( sum( ( data{i}.y - kooplin{i}.y ).^2 , 2 ) );
    kooplin{i}.err_ave = sum( kooplin{i}.err ) / length(kooplin{i}.err);
    kooplin{i}.err_cm = in2cm( kooplin{i}.err );
    kooplin{i}.err_ave_cm = in2cm( kooplin{i}.err_ave );
end

%% simulate linear state space model

% load in the linear state space model
load('larmv5_linSS.mat');
lin = cell( size(data) );
for i = 1 : length(data)
    lin{i}.t = data{i}.t;
    lin{i}.u = data{i}.u;
    
%     % scale down the data
%     datasc.t = data{i}.t;
%     datasc.u = ksysid_lin.scaledown.u( data{i}.u );
%     datasc.y = ksysid_lin.scaledown.y( data{i}.y );
    
    lin{i}.x = zeros( length( data{i}.t ) , 4 );    % preallocate
    lin{i}.x(1,:) = [ data{i}.y(1,1) , data{i}.y(1,1) , data{i}.y(1,2) , data{i}.y(1,2) ];  % initial state
    lin{i}.y = zeros( length( data{i}.t ) , 2 );    % preallocate
    lin{i}.y(1,:) = data{i}.y(1,:); % initial output
    
    % simulate the system
    for j = 1 : length( kooplin{i}.t ) - 1
        lin{i}.x(j+1,:) = ( model.A * lin{i}.x(j,:)' + model.B * data{i}.u(j,:)' )';
        lin{i}.y(j+1,:) = ( model.C * lin{i}.x(j+1,:)' )';
    end
    
%     % scale data back up
%     kooplin{i}.u = data{i}.u;
%     kooplin{i}.y = ksysid_lin.scaleup.y( kooplin{i}.y );
    
    % distance error at each timestep
    lin{i}.err = sqrt( sum( ( data{i}.y - lin{i}.y ).^2 , 2 ) );
    lin{i}.err_ave = sum( lin{i}.err ) / length(lin{i}.err);
    lin{i}.err_cm = in2cm( lin{i}.err );
    lin{i}.err_ave_cm = in2cm( lin{i}.err_ave );
end

%% simulate nonlinear Koopman model

% % load in the nonlinear koopman model (continuous)
% matcontents = load('n-2_m-3_del-0_2020-04-29_21-20.mat');
% ksysid = matcontents.sysid_class;
% koop = cell( size(data) );
% for i = 1 : length(data)
% %     Ts_nl = ( data{i}.t(end) + 1 ) / length( data{i}.t );   % actual average timestep
%     koop{i}.t = data{i}.t;
%     koop{i}.u = ksysid.scaledown.u( data{i}.u );
%     ykoop0 = ksysid.scaledown.y( data{i}.y(1,:) );
%     [ ~ , koop{i}.y ] = ode45( @(t,y) ksysid.model.F_func( y , koop{i}.u( min( floor(t/Ts)+1 , size(koop{i}.u,1) ) ,:)' ) , koop{i}.t , ykoop0' );
%     
%     % scale data back up
%     koop{i}.u = data{i}.u;
%     koop{i}.y = ksysid.scaleup.y( koop{i}.y );
%     
%     % distance error at each timestep
%     koop{i}.err = sqrt( sum( ( data{i}.y - koop{i}.y ).^2 , 2 ) );
%     koop{i}.err_ave = sum( koop{i}.err ) / length(koop{i}.err);
%     koop{i}.err_cm = in2cm( koop{i}.err );
%     koop{i}.err_ave_cm = in2cm( koop{i}.err_ave );
% end

% load in the nonlinear koopman model (discrete)
matcontents = load('nonlinear-discrete-laserrobot-model.mat');
ksysid = matcontents.temp;
koop = cell( size(data) );
for i = 1 : length(data)
    koop{i}.t = data{i}.t;
    
    % scale down the data (re-uses scaling function from linear model)
    datasc.t = data{i}.t;
    datasc.u = ksysid_lin.scaledown.u( data{i}.u );
    datasc.y = ksysid_lin.scaledown.y( data{i}.y );
    
    koop{i}.zeta = zeros( length(koop{i}.t) , size( kooplin{1}.zeta , 2 ) );   % preallocate  
    koop{i}.zeta(1,:) = kooplin{i}.zeta(1,:);   % initial value
    koop{i}.y = zeros( size( data{i}.y ) );  % preallocate
    koop{i}.y(1,:) = datasc.y(1,:);  % initial value
    
    % simulate the (scaled) system
    for j = 1 : length( koop{i}.t ) - 1
        koop{i}.zeta(j+1,:) = ksysid.F_func( koop{i}.zeta(j,:)' , datasc.u(j,:)' )';
        koop{i}.y(j+1,:) = koop{i}.zeta(j+1,1:2);
    end
    
    % scale data back up
    koop{i}.u = data{i}.u;
    koop{i}.y = ksysid_lin.scaleup.y( koop{i}.y );
    
    % distance error at each timestep
    koop{i}.err = sqrt( sum( ( data{i}.y - koop{i}.y ).^2 , 2 ) );
    koop{i}.err_ave = sum( koop{i}.err ) / length(koop{i}.err);
    koop{i}.err_cm = in2cm( koop{i}.err );
    koop{i}.err_ave_cm = in2cm( koop{i}.err_ave );
end

%% simulate neural network model

nnet = cell( size(data) );
for i = 1 : length(data)
    nnet{i}.t = data{i}.t;
    nnet{i}.u = data{i}.u;

    nnet{i}.y = sim_nnet_model( data{i}.u , data{i}.y(1,:) , data{i}.u(1,:) , data{i}.y(1,:) );

%     [ ynnet , ~ ] = nnet_model_nio( nnet{i}.u' , [0 0 0]' );
%     nnet{i}.y = ynnet'; 
    
    % distance error at each timestep
    nnet{i}.err = sqrt( sum( ( data{i}.y - nnet{i}.y ).^2 , 2 ) );
    nnet{i}.err_ave = sum( nnet{i}.err ) / length(nnet{i}.err);
    nnet{i}.err_cm = in2cm( nnet{i}.err );
    nnet{i}.err_ave_cm = in2cm( nnet{i}.err_ave );
end


%% simulate NLARX model

% load in the NLARX model from file
load('nlarx1.mat');
% load('nlarx7.mat');
% nlarx1 = nlarx7;
nlarx = cell( size(data) );
for i = 1 : length( data )
    nlarx{i}.t = data{i}.t;
    nlarx{i}.u = data{i}.u;
    udata = iddata( [] , nlarx{i}.u , nlarx1.Ts ); 
    ydata = sim( nlarx1 , udata );
    nlarx{i}.y = ydata.y;
    
%     nlarx{i}.y(1,:) = data{i}.y(1,:);
%     for j = 1 : length( data{i}.t ) - 1
%         IC.Input = nlarx{i}.u(1:j,:);
%         IC.Output = nlarx{i}.y(1:j,:);
%         opt = simOptions('InitialCondition',IC);
%         nlarx{i}.y(j+1,:) = sim( nlarx1 , nlarx{i}.u(j,:) );
%     end
    
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

%% create the results output

res = struct;
res.real = data;    % actual system
res.koop = koop;    % nonlinear koopman
res.nnet = nnet;    % neural network
res.nlarx = nlarx; 
res.nlhw = nlhw;
res.kooplin = kooplin;  % linear koopman model
res.lin = lin;          % linear state space model (not lifted)

%% evaluate the total average error for each model

% sum the error from all the trials
koop_allerr = [];
nnet_allerr = [];
nlarx_allerr = [];
nlhw_allerr = [];
kooplin_allerr = [];
lin_allerr = [];

for i = 1 : length( data )
    % concatenate all the error vectors
    koop_allerr = [ koop_allerr ; koop{i}.err ];
    nnet_allerr = [ nnet_allerr ; nnet{i}.err ];
    nlarx_allerr = [ nlarx_allerr ; nlarx{i}.err ];
    nlhw_allerr = [ nlhw_allerr ; nlhw{i}.err ];
    kooplin_allerr = [ kooplin_allerr ; kooplin{i}.err ];
    lin_allerr = [ lin_allerr ; lin{i}.err ];
end

% calculate average error for each model
totsteps = size( koop_allerr , 1 );
err.koop.ave = sum( koop_allerr ) / totsteps;
err.nnet.ave = sum( nnet_allerr ) / totsteps;
err.nlarx.ave = sum( nlarx_allerr ) / totsteps;
err.nlhw.ave = sum( nlhw_allerr ) / totsteps;
err.kooplin.ave = sum( kooplin_allerr ) / totsteps;
err.lin.ave = sum( lin_allerr ) / totsteps;

% calculate standard deviation
err.koop.std = std( koop_allerr );
err.nnet.std = std( nnet_allerr );
err.nlarx.std = std( nlarx_allerr );
err.nlhw.std = std( nlhw_allerr );
err.kooplin.std = std( kooplin_allerr );
err.lin.std = std( lin_allerr );



%% plot the error comparison over all validation trials

% Define custom colors:
cb(1,:) = [69,117,180]/255;
cb(2,:) = [145,191,219]/255;
cb(3,:) = [254,224,144]/255;
cb(4,:) = [252,141,89]/255;
cb(5,:) = [215,48,39]/255;

% Define more colors
cb_red = [228,26,28] ./ 255;   % red (for linear state space)
cb_blue = [55,126,184] ./ 255;  % blue (for linear koopman)
cb_orange = [255,127,0] ./ 255;   % orange (for nonlinear koopman)
cb_grey = [124,124,124] ./ 255;
cb_brown = [166,86,40] ./ 255;

% bar chart showing TOTAL NRMSE across all states and all trials
figure
hold on
bar(1 , err.kooplin.ave * 2.54 , 'FaceColor' , cb_blue );
bar(2 , err.lin.ave * 2.54 , 'FaceColor' , cb_red );
bar(3 , err.koop.ave * 2.54 , 'FaceColor' , cb_orange );
bar(4 , err.nnet.ave * 2.54 , 'FaceColor' , cb_brown );
bar(5 , err.nlarx.ave * 2.54 , 'FaceColor' , cb_grey );
bar(6 , err.nlhw.ave * 2.54 , 'FaceColor' , cb_grey );
xticks(1:6);
xtickangle(45);
xticklabels( {'Koopman (Linear)', 'State Space (Linear)' , 'Koopman (Nonlinear)' , 'Neural Network', 'NLARX', 'Ham.-Wiener'} );
ylabel('Average Error (cm)');
ylim([0 12]);
% eb = errorbar(1:numSystems-1, [mean_NRMSE_alltrials(1:3)' * 100, mean_NRMSE_alltrials(5:end)' * 100], [std_NRMSE_alltrials(1:3)' * 100, std_NRMSE_alltrials(5:end)' * 100], '.', 'CapSize', 14, 'LineWidth', 1.5, 'Color', 'k');
eb = errorbar(1:6, [err.kooplin.ave , err.lin.ave , err.koop.ave , err.nnet.ave , err.nlarx.ave , err.nlhw.ave] * 2.54, [err.kooplin.std , err.lin.std , err.koop.std , err.nnet.std , err.nlarx.std , err.nlhw.std] * 2.54, '.', 'CapSize', 14, 'LineWidth', 1.5, 'Color', 'k');
box on
hold off


%% plot the error comparison over all validation trials (excluding NLARX and Hamm-Weiner)

% Define custom colors:
cb(1,:) = [69,117,180]/255;
cb(2,:) = [145,191,219]/255;
cb(3,:) = [254,224,144]/255;
cb(4,:) = [252,141,89]/255;
cb(5,:) = [215,48,39]/255;

% Define more colors
cb_red = [228,26,28] ./ 255;   % red (for linear state space)
cb_blue = [55,126,184] ./ 255;  % blue (for linear koopman)
cb_orange = [255,127,0] ./ 255;   % orange (for nonlinear koopman)
cb_grey = [150,150,150] ./ 255;
cb_brown = [166,86,40] ./ 255;

% bar chart showing TOTAL NRMSE across all states and all trials
figure
hold on
bar(1 , err.kooplin.ave * 2.54 , 'FaceColor' , cb_blue );
bar(2 , err.lin.ave * 2.54 , 'FaceColor' , cb_red );
bar(4 , err.koop.ave * 2.54 , 'FaceColor' , cb_orange );
bar(5 , err.nnet.ave * 2.54 , 'FaceColor' , cb_brown );
xticks([0.85 1.15 2 3.85 4.15 5]);
xtickangle(45);
xticklabels( {'Koopman' , '(Linear)     ' , 'State Space' , 'Koopman' , '(Nonlinear)   ' , 'Neural Network'} );
ylabel('Average Error (cm)');
ylim([0 8]);
% eb = errorbar(1:numSystems-1, [mean_NRMSE_alltrials(1:3)' * 100, mean_NRMSE_alltrials(5:end)' * 100], [std_NRMSE_alltrials(1:3)' * 100, std_NRMSE_alltrials(5:end)' * 100], '.', 'CapSize', 14, 'LineWidth', 1.5, 'Color', 'k');
eb = errorbar([1 2 4 5], [err.kooplin.ave , err.lin.ave , err.koop.ave , err.nnet.ave] * 2.54, [err.kooplin.std , err.lin.std , err.koop.std , err.nnet.std] * 2.54, '.', 'CapSize', 14, 'LineWidth', 1.5, 'Color', 'k');
box on;
set(gca, 'YGrid', 'on', 'XGrid', 'off');
hold off

%% plot the error comparison over all validation trials (excluding NLARX and Hamm-Weiner) but using colors for dissertation defense

% Define color map
colormap lines;
cb = colormap;

% % Define custom colors:
% cb(1,:) = [69,117,180]/255;
% cb(2,:) = [145,191,219]/255;
% cb(3,:) = [254,224,144]/255;
% cb(4,:) = [252,141,89]/255;
% cb(5,:) = [215,48,39]/255;

% Define more colors
cb_blue = cb(1,:);  % blue (for linear koopman)
cb_orange = cb(2,:);   % orange (for nonlinear koopman)


% bar chart showing TOTAL NRMSE across all states and all trials
figure
hold on
bar(1 , err.kooplin.ave * 2.54 , 'FaceColor' , cb_blue );
bar(2 , err.lin.ave * 2.54 , 'FaceColor' , cb_grey );
bar(4 , err.koop.ave * 2.54 , 'FaceColor' , cb_orange );
bar(5 , err.nnet.ave * 2.54 , 'FaceColor' , cb_grey );
xticks([0.85 1.15 2 3.85 4.15 5]);
xtickangle(45);
xticklabels( {'Koopman' , '(Linear)     ' , 'State Space' , 'Koopman' , '(Nonlinear)   ' , 'Neural Network'} );
ylabel('Average Error (cm)');
ylim([0 5]);
% eb = errorbar(1:numSystems-1, [mean_NRMSE_alltrials(1:3)' * 100, mean_NRMSE_alltrials(5:end)' * 100], [std_NRMSE_alltrials(1:3)' * 100, std_NRMSE_alltrials(5:end)' * 100], '.', 'CapSize', 14, 'LineWidth', 1.5, 'Color', 'k');
eb = errorbar([1 2 4 5], [err.kooplin.ave , err.lin.ave , err.koop.ave , err.nnet.ave] * 2.54, [err.kooplin.std , err.lin.std , err.koop.std , err.nnet.std] * 2.54, '.', 'CapSize', 14, 'LineWidth', 1.5, 'Color', 'k');
box on;
set(gca, 'YGrid', 'on', 'XGrid', 'off');
hold off



