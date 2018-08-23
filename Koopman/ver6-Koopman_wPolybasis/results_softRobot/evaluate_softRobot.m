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

numChops = 3;   % slices to chop data into get 10s trials

for i = 1 : params.numVals
    valID = ['val', num2str(i)];
    [t(i), x(i), u(i)] = chop_data(numChops, data.(valID)); % chop data into 10s chunks
end


%% Simulate Koopman model over all validation sets

valCount = 0;   % counter of total validation trials
for j = 1 : params.numVals
    for k = 1 : numChops
        valCount = valCount + 1;    % increment validation trial count       
        
        % simulate koopman system
        chopID = ['s', num2str(k)];
        tval = t(i).(chopID);
        xval = x(i).(chopID);
        uval = u(i).(chopID);
        x0sim = x(i).(chopID)(1,:)'; % same initial state as validation data initial state
        [tsysid, xsysid] = ode45(@(t,x) vf_koopman(x, get_u(t, tval, uval)), tval, x0sim);
    
        % save results as iddata object
        valID = ['z', num2str(valCount)];
        val.(valID) = iddata( xval, uval, Ts, 'Name', 'Actual');      % actual system
        koop.(valID) = iddata( xsysid, uval, Ts, 'Name', 'Koopman');   % koopman system
    end
end

%% Simulate Neural Network model over all validation sets

% STILL NEED TO DO THIS PART

%% Compare performance of all systems

load('MatlabSysid.mat');    % loads in models from Matlab sysid toobox
fit = cell(valCount, 1);
for i = 1 : valCount
    valID = ['z', num2str(i)];  % trial identifier
%     x0 = val.(valID).y(1,:)';
%     opt = compareOptions('InitialCondition', x0);
    [yh, fit{i}, x0] = compare(val.(valID), koop.(valID), ARMAX, HammersteinWiener, NLARX, StateSpace);
end

% FIGURE OUT A WAY OF AVERAGING ERROR FROM ALL THESE SYSTEMS!!

%% THIS IS WHERE I LEFT OFF




error = struct;     % error results comparing real and koopman system
koopsim = struct;   % simulation results for koopman system

for j = 1 : valparams.numVals
    %% simulate the behavior of learned system
    
    % isolate the jth validation trial
    valID = ['val', num2str(j)];
    valdata = data.(valID);
    
    tspan = valdata.t;
    
    x0sim = valdata.x(1,:)'; % same initial state as validation data initial state
    [tsysid, xsysid] = ode45(@(t,x) vf_koopman(x, get_u(t, x, valdata, valparams)), tspan, x0sim);
    
    % simulated forward using the transpose of Koopman operator
    xselector = [zeros(valparams.n,1), eye(valparams.n), zeros(valparams.n, valparams.N - valparams.n - 1)]; % matrix to extract state from lifted state
    xkoop = zeros(length(tspan), valparams.n);
    xkoop(1,:) = x0sim';
    for i = 2 : length(tspan)
        ti = tspan(i);
        xnext = xselector * koopman.U' * polyLift( xkoop(i-1,:)' , get_u(ti, 0, valdata, valparams) );
        xkoop(i,:) = xnext';
    end
    
    % real system
    [treal, xreal] = deal(tspan, valdata.x(1:length(tspan) , :));
    
    
    %% quantify the error between real behavior and simulated behavior
    
    terror = treal;
    xerror = abs( xreal - xsysid );
    xerrormax = max(max(xerror(:,1:ceil(valparams.n/2))));
    % xerrormin = min(min(xerror(:,1:ceil(valparams.n/2))));
    RMSE = sqrt( sum( (xreal - xsysid).^2 ) / length(terror) );
    
    % defind output
    error.(valID).terror = terror;
    error.(valID).xerror = xerror;
    error.(valID).RMSE = RMSE;
    
    koopsim.(valID).t = tsysid;
    koopsim.(valID).x = xsysid;
    
    
    %% plot the results
    
    if valparams.ploton
        figure
        subplot(4,1,1)
        plot(treal, xreal(:,1:ceil(valparams.n/2)))
        title('Real system')
        subplot(4,1,2)
        plot(tsysid, xsysid(:,1:ceil(valparams.n/2)))
        title('Identified system')
        subplot(4,1,3)
        hold on
        plot(terror, xerror(:,1:ceil(valparams.n/2)))
        plot(terror, xerrormax * ones(size(terror)), '--')
        title('Error')
        hold off
        subplot(4,1,4)
        plot(tspan, xkoop(:,1:ceil(valparams.n/2)))
        title('Koopman Transpose Sim')
    end
    
    
    % % animate the results
    % animate_doublePendulum(sol_real, sol_sysid, valparams);

end


end


function u = get_u(t, tval, uval)
%get_u: Interpolates to estimate the value of the input at a specific t

u = interp1(tval, uval, t)';

end