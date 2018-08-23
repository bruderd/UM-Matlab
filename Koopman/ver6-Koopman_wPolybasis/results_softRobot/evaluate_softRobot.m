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
        [xnnet ,~] = NLinout_NeuralNetworkFunction(uval', zeros(3,2));
        
        % save results as iddata object
        nnet.(valID) = iddata( xnnet', uval, Ts, 'Name', 'NeuralNetwork' );   % neural network system
        
    end
end

%% Simulate Neural Network model over all validation sets

% STILL NEED TO DO THIS PART


%% Compare performance of all systems

% load('MatlabSysid.mat');    % loads in models from Matlab sysid toobox
fit = cell(valCount, 1);
for i = 1 : valCount
%     sysID = matSystems{i}.name;
    valID = ['z', num2str(i)];  % trial identifier
    x0 = val.(valID).y(1,:)';
    opt = compareOptions('InitialCondition', 'z');
%     [yh, fit{i}, x0] = compare(val.(valID), koop.(valID), matSys{1}.(valID), matSys{2}.(valID), matSys{3}.(valID), matSys{4}.(valID));
    [yh, fit{i}, x0] = compare(val.(valID), koop.(valID), nnet.(valID), stateSpace.(valID), matSystems{1}, matSystems{2}, matSystems{3}, opt);
end

% Compute average %FIT across all states and all trials for each system
for i = 1 : numSystems
    for j = 1 : valCount
        states_fit(j,:) = fit{j}{i}(:)';    % NRMSE for all of the states 
%         mean_states_fit(j,i) = sum( fit{j}{i}(:) ) / params.n;  % mean of error across all states in jth trial
    end
    mean_states_fit(i,:) = mean( states_fit, 1 );
    std_states_fit(i,:) = std( states_fit, 0, 1 );
    
    mean_fit(i) = mean2( states_fit );
    std_fit(i) = std2( states_fit );
%     mean_fit(i) = sum( mean_states_fit(:,i) ) / valCount;   % mean of error across all states and trials
%     std_fit(i) = std( mean_states_fit(:,i) );   % standard deviation accross all states and trials
end


%% Plot the results

% figure
% bar([mean_states_fit(1:3,:); mean_states_fit(5:end,:)] );
% errorbar(1:5, [mean_states_fit(1:3,:); mean_states_fit(5:end,:)], [std_states_fit(1:3,:); std_states_fit(5:end,:)]);

figure
hold on
b = bar([mean_fit(1:3), mean_fit(5:end)], 'b');
xticks(1:5);
xtickangle(60);
xticklabels( {'Koopman', 'Neural Network', 'State Space', matSystems{2}.name, matSystems{3}.name} );
ylabel('Fit (%)');
eb = errorbar(1:numSystems-1, [mean_fit(1:3), mean_fit(5:end)], [std_fit(1:3), std_fit(5:end)], '.', 'CapSize', 14, 'LineWidth', 2, 'Color', 'k');
hold off




end


function u = get_u(t, tval, uval)
%get_u: Interpolates to estimate the value of the input at a specific t

u = interp1(tval, uval, t)';

end