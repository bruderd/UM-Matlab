function [error, koopsim] = evaluate_softRobot( data )
%evaluate_softRobot: Evaluates the performance of koopman model for soft
%robot system verses Matlab sysid toolbox models and neural network
%   Detailed explanation goes here

%% parameter values that need to be set


%% If data argument not supplied, prompt uset to load data from file
if ~exist('data', 'var')
    % Prompt user to identify data file
    [data_file,data_path] = uigetfile('..\dataFiles\*.mat');
    matcontents = load([data_path, data_file]); % must be a .mat file
    data = matcontents.data; 
end

% extract some values
params = data.valparams;
Ts = params.Ts;

%% Split data into 10s validation steps

numChops = 1;   % slices to chop data into get 10s trials

for i = 1 : params.numVals
    valID = ['val', num2str(i)];
    [t(i), x(i), u(i)] = chop_data(numChops, data.(valID)); % chop data into 10s chunks
end


%% Simulate Koopman model over all validation sets

S = load('MatlabSysid.mat');    % loads in models from Matlab sysid toobox
matSystems = struct2cell(S);    % converts struct into indexable thing
numSystems = size(matSystems,1) + 2; % number of systems to be compared

% % convert any polynomial models to state space models
% for i = 1 : numSystems - 1
%     if strcmp( class(matSystems{i}), 'idpoly')
%         matSystems{i} = ss(matSystems{i});  % convert to ss model so ic can be set
%     end
% end

valCount = 0;   % counter of total validation trials
for j = 1 : params.numVals
    for k = 1 : numChops
        valCount = valCount + 1;    % increment validation trial count       
        
        % simulate koopman system
        chopID = ['s', num2str(k)];
        tval = t(j).(chopID);
        xval = x(j).(chopID);
        uval = u(j).(chopID);
        x0sim = x(j).(chopID)(1,:)'; % same initial state as validation data initial state
        [tkoop, xkoop] = ode45(@(t,x) vf_koopman(x, get_u(t, tval, uval)), tval, x0sim);
    
        % save results as iddata object
        valID = ['z', num2str(valCount)];
        time.(valID) = tval;   % save the time vector for this trial
        val.(valID) = iddata( xval, uval, Ts, 'Name', 'Actual');      % actual system
        koop.(valID) = iddata( xkoop, uval, Ts, 'Name', 'Koopman');   % koopman system
        
        % simulate state space model using real IC
        for i = 1 : numSystems - 2
            if strcmp( class(matSystems{i}) , 'idss')
                opt = simOptions('InitialCondition', 'z');
                xss = sim(matSystems{i}, uval, opt);
                
%               % shift result to make initial conditions line up
%               xmatsys = xmatsys + ( x0sim' - xmatsys(1,:) );
                
                % save results as iddata object
                sysID = matSystems{i}.name;
                stateSpace.(valID) = iddata( xss, uval, Ts, 'Name', sysID );   % state space system
            end
        end
        
        % simulate neural network model using real IC
        [xnnet ,~] = NLinout_NeuralNetworkFunction_v6(uval', zeros(3,2));
        
        % save results as iddata object
        nnet.(valID) = iddata( xnnet', uval, Ts, 'Name', 'NeuralNetwork' );   % neural network system
        
    end
end


%% Compare performance of all systems

% load('MatlabSysid.mat');    % loads in models from Matlab sysid toobox
fit = cell(valCount, 1);
yh = cell(valCount, 1);
for i = 1 : valCount
%     sysID = matSystems{i}.name;
    valID = ['z', num2str(i)];  % trial identifier
    x0 = val.(valID).y(1,:)';
    opt = compareOptions('InitialCondition', 'z');
%     [yh, fit{i}, x0] = compare(val.(valID), koop.(valID), matSys{1}.(valID), matSys{2}.(valID), matSys{3}.(valID), matSys{4}.(valID));
    [yh{i}, fit{i}, x0] = compare(val.(valID), koop.(valID), nnet.(valID), stateSpace.(valID), matSystems{1}, matSystems{2}, matSystems{3}, opt);
end

% Compute average FIT across all states and all trials for each system
for i = 1 : numSystems
    for j = 1 : valCount
        states_fit(j,:,i) = fit{j}{i}(:)';    % NRMSE for all of the states 
%         mean_states_fit(j,i) = sum( fit{j}{i}(:) ) / params.n;  % mean of error across all states in jth trial
    end
    mean_states_fit(i,:) = mean( states_fit(:,:,i), 1 );
    std_states_fit(i,:) = std( states_fit(:,:,i), 0, 1 );
    
    mean_fit(i) = mean2( states_fit(:,:,i) );
    std_fit(i) = std2( states_fit(:,:,i) );
%     mean_fit(i) = sum( mean_states_fit(:,i) ) / valCount;   % mean of error across all states and trials
%     std_fit(i) = std( mean_states_fit(:,i) );   % standard deviation accross all states and trials
end


% Compute NRMSE (Not 1-NRMSE like matlab does) for each trial
for i = 1 : numSystems
    for j = 1 : valCount
        valID = ['z', num2str(j)];  % trial identifier
        states_NRMSE(j,:,i) = sqrt( sum( ( yh{j}{i}.y  - val.(valID).y ).^2 ) / size(val.(valID).y, 1) ) ./ ( max(data.alltrials.x, [], 1) - min(data.alltrials.x, [], 1) ); % NRMSE (normalized by max value of each state)
%         states_NRMSE(j,:,i) = sqrt( sum( ( yh{j}{i}.y  - val.(valID).y ).^2 ) ./ size(val.(valID).y, 1) ) ./ abs( mean( val.(valID).y, 1 ) ); % NRMSE (normalized by absolute value of trial mean), bad becasuse divides by almost zero!!
    end
    mean_states_NRMSE(i,:) = mean( states_NRMSE(:,:,i), 1 );
    std_states_NRMSE(i,:) = std( states_NRMSE(:,:,i), 0, 1 );
    
    mean_NRMSE(i) = mean2( states_NRMSE(:,:,i) );
    std_NRMSE(i) = std2( states_NRMSE(:,:,i) );
end

% Compute RMSE (in meters) for each trial
for i = 1 : numSystems
    for j = 1 : valCount
        valID = ['z', num2str(j)];  % trial identifier
%         states_RMSE(j,:,i) = -100*( ( fit{j}{i}(:)' ./ 100 ) - 1 ) * norm( val.(valID).y - mean( val.(valID).y, 1 ) ) ./size(val.(valID).y, 1);    % RMSE for all of the states (BAD)
        states_RMSE(j,:,i) = sqrt( sum( ( yh{j}{i}.y  - val.(valID).y ).^2 ) ./ size(val.(valID).y, 1) );  % RMSE for all of the states 
    end
    mean_states_RMSE(i,:) = mean( states_RMSE(:,:,i), 1 );
    std_states_RMSE(i,:) = std( states_RMSE(:,:,i), 0, 1 );
    
    mean_RMSE(i) = mean2( states_RMSE(:,:,i) );
    std_RMSE(i) = std2( states_RMSE(:,:,i) );
end

% Compute total RMSE and NRMSE across all trials!
states_squerror = zeros(numSystems, params.n);
bigN = 0;
for i = 1 : numSystems
    for j = 1 : valCount
        valID = ['z', num2str(j)];  % trial identifier
        states_squerror(i,:) = states_squerror(i,:) + sum( ( yh{j}{i}.y  - val.(valID).y ).^2 ) ;  % adding the square error from all trials
        bigN = bigN + size(val.(valID).y, 1) ;  % total number of data points for all trials
    end
end
states_RMSE_alltrials = sqrt( states_squerror ./ bigN );
states_NRMSE_alltrials = states_RMSE_alltrials ./ ( max(data.alltrials.x, [], 1) - min(data.alltrials.x, [], 1) );
mean_RMSE_alltrials = mean(states_RMSE_alltrials, 2);
mean_NRMSE_alltrials = mean(states_NRMSE_alltrials, 2);
std_RMSE_alltrials = std(states_RMSE_alltrials, [], 2);
std_NRMSE_alltrials = std(states_NRMSE_alltrials, [], 2);


%% Plot the results

% Define custom colors:
cb(1,:) = [69,117,180]/255;
cb(2,:) = [145,191,219]/255;
cb(3,:) = [254,224,144]/255;
cb(4,:) = [252,141,89]/255;
cb(5,:) = [215,48,39]/255;


% bar chart showing average fit across all states and all trials
figure
hold on
b = bar([mean_fit(1:3), mean_fit(5:end)], 'b');
xticks(1:5);
xtickangle(60);
xticklabels( {'Koopman', 'Neural Network', 'State Space', matSystems{2}.name, matSystems{3}.name} );
ylabel('1 - NRMSE (%)');
eb = errorbar(1:numSystems-1, [mean_fit(1:3), mean_fit(5:end)], [std_fit(1:3), std_fit(5:end)], '.', 'CapSize', 14, 'LineWidth', 2, 'Color', 'k');
hold off

% bar chart showing average NRMSE across all states and all trials
figure
hold on
b = bar([mean_NRMSE(1:3), mean_NRMSE(5:end)], 'b');
xticks(1:5);
xtickangle(60);
xticklabels( {'Koopman', 'Neural Network', 'State Space', matSystems{2}.name, matSystems{3}.name} );
ylabel('NRMSE (%)');
eb = errorbar(1:numSystems-1, [mean_NRMSE(1:3), mean_NRMSE(5:end)], [std_NRMSE(1:3), std_NRMSE(5:end)], '.', 'CapSize', 14, 'LineWidth', 2, 'Color', 'k');
hold off

% bar chart showing TOTAL NRMSE across all states and all trials
figure
hold on
for i = 1 : 3
    bar(i, mean_NRMSE_alltrials(i) * 100, 'FaceColor', cb(i,:));
end
for i = 5 : numSystems
    bar(i-1, mean_NRMSE_alltrials(i) * 100, 'FaceColor', cb(i-1,:));
end
xticks(1:5);
xtickangle(60);
xticklabels( {'Koopman', 'Neural Network', 'State Space', 'Hammerstein-Wiener', matSystems{3}.name} );
ylabel('NRMSE (%)');
ylim([0 10]);
eb = errorbar(1:numSystems-1, [mean_NRMSE_alltrials(1:3)' * 100, mean_NRMSE_alltrials(5:end)' * 100], [std_NRMSE_alltrials(1:3)' * 100, std_NRMSE_alltrials(5:end)' * 100], '.', 'CapSize', 14, 'LineWidth', 1.5, 'Color', 'k');
box on
hold off

% bar chart showing average RMSE across all states and all trials
figure
hold on
b = bar([mean_RMSE(1:3), mean_RMSE(5:end)], 'b');
xticks(1:5);
xtickangle(60);
xticklabels( {'Koopman', 'Neural Network', 'State Space', matSystems{2}.name, matSystems{3}.name} );
ylabel('RMSE (m)');
eb = errorbar(1:numSystems-1, [mean_RMSE(1:3), mean_RMSE(5:end)], [std_RMSE(1:3), std_RMSE(5:end)], '.', 'CapSize', 14, 'LineWidth', 2, 'Color', 'k');
hold off

% bar chart showing average fit across all trials
figure
hold on
b2 = bar([mean(states_fit(:,:,1),1); mean(states_fit(:,:,2),1); mean(states_fit(:,:,3),1); mean(states_fit(:,:,5),1); mean(states_fit(:,:,6),1)]);
xticks(1:5);
xtickangle(60);
xticklabels( {'Koopman', 'Neural Network', 'State Space', matSystems{2}.name, matSystems{3}.name} );
ylabel('1 - NRMSE (%)');
% eb2 = errorbar(1:numSystems-1, [mean_fit(1:3), mean_fit(5:end)], [std_fit(1:3), std_fit(5:end)], '.', 'CapSize', 14, 'LineWidth', 2, 'Color', 'k');
hold off


% bar chart showing average NRMSE across all trials
figure
hold on
b4 = bar([mean(states_NRMSE(:,:,1),1); mean(states_NRMSE(:,:,2),1); mean(states_NRMSE(:,:,3),1); mean(states_NRMSE(:,:,5),1); mean(states_NRMSE(:,:,6),1)]);
xticks(1:5);
xtickangle(60);
xticklabels( {'Koopman', 'Neural Network', 'State Space', matSystems{2}.name, matSystems{3}.name} );
ylabel('NRMSE (%)');
% plot([0, 6], [1, 1], 'k--')
% eb4 = errorbar(1:numSystems-1, [mean_fit(1:3), mean_fit(5:end)], [std_fit(1:3), std_fit(5:end)], '.', 'CapSize', 14, 'LineWidth', 2, 'Color', 'k');
hold off

% bar chart showing average RMSE across all trials
figure
hold on
b5 = bar([mean(states_RMSE(:,:,1),1); mean(states_RMSE(:,:,2),1); mean(states_RMSE(:,:,3),1); mean(states_RMSE(:,:,5),1); mean(states_RMSE(:,:,6),1)]);
xticks(1:5);
xtickangle(60);
xticklabels( {'Koopman', 'Neural Network', 'State Space', matSystems{2}.name, matSystems{3}.name} );
ylabel('RMSE (m)');
% eb4 = errorbar(1:numSystems-1, [mean_fit(1:3), mean_fit(5:end)], [std_fit(1:3), std_fit(5:end)], '.', 'CapSize', 14, 'LineWidth', 2, 'Color', 'k');
hold off


%% Plot XYZ trajectories of all the systems for each validation trial

for i = 1 : valCount
    valID = ['z', num2str(i)];  % trial identifier
    
    figure;
    subplot(3,1,1)
        hold on
        plot(time.(valID), val.(valID).y(:,1) * 100, ':k', 'LineWidth', 2);
        plot(time.(valID), yh{i}{1}.y(:,1) * 100, 'Color', cb(1,:), 'LineWidth', 2);
%         plot(time.(valID), yh{i}{2}.y(:,1) * 100, 'Color', cb(2,:), 'LineWidth', 3);
%         plot(time.(valID), yh{i}{3}.y(:,1) * 100, 'Color', cb(3,:), 'LineWidth', 3);
%         plot(time.(valID), yh{i}{5}.y(:,1) * 100, 'Color', cb(4,:), 'LineWidth', 3);
%         plot(time.(valID), yh{i}{6}.y(:,1) * 100, 'Color', cb(5,:), 'LineWidth', 3);
%         ylabel('$x_1$ (cm)', 'Interpreter', 'Latex')
        ylabel('x_1')
        xticklabels([])
        hold off
        box on
    subplot(3,1,2)
        hold on
        plot(time.(valID), val.(valID).y(:,2) * 100, ':k', 'LineWidth', 2);
        plot(time.(valID), yh{i}{1}.y(:,2) * 100, 'Color', cb(1,:), 'LineWidth', 2);
%         plot(time.(valID), yh{i}{2}.y(:,2) * 100, 'Color', cb(2,:), 'LineWidth', 3);
%         plot(time.(valID), yh{i}{3}.y(:,2) * 100, 'Color', cb(3,:), 'LineWidth', 3);
%         plot(time.(valID), yh{i}{5}.y(:,2) * 100, 'Color', cb(4,:), 'LineWidth', 3);
%         plot(time.(valID), yh{i}{6}.y(:,2) * 100, 'Color', cb(5,:), 'LineWidth', 3);
%         ylabel('$x_2$ (cm)', 'Interpreter', 'Latex')
        ylabel('x_2')
        xticklabels([])
        hold off
        box on
    subplot(3,1,3)
        hold on
        plot(time.(valID), val.(valID).y(:,3) * 100, ':k', 'LineWidth', 2);
        plot(time.(valID), yh{i}{1}.y(:,3) * 100, 'Color', cb(1,:), 'LineWidth', 2);
%         plot(time.(valID), yh{i}{2}.y(:,3) * 100, 'Color', cb(2,:), 'LineWidth', 3);
%         plot(time.(valID), yh{i}{3}.y(:,3) * 100, 'Color', cb(3,:), 'LineWidth', 3);
%         plot(time.(valID), yh{i}{5}.y(:,3) * 100, 'Color', cb(4,:), 'LineWidth', 3);
%         plot(time.(valID), yh{i}{6}.y(:,3) * 100, 'Color', cb(5,:), 'LineWidth', 3);
%         ylabel('$x_3$ (cm)', 'Interpreter', 'Latex')
        ylabel('x_3')
        xlabel('Time (s)')
        hold off
        box on
        legend({'Real System', 'Koopman Model'}, 'Location', 'NorthEast')
end

end


function u = get_u(t, tval, uval)
%get_u: Interpolates to estimate the value of the input at a specific t

u = interp1(tval, uval, t)';

end