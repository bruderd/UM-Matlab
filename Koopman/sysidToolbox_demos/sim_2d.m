%% Simulate 2D system
clear;


tspan = [0 50];
x0 = [0, 0];
options = odeset('OutputFcn',@(t,x,flag) input_2d(t, x, flag), 'refine', 1);

input_2d(0,0,'clear');
[t, y] = ode45(@(t,x) myode_1d(t, x, input_2d(t, x, 'pass2ode')), tspan, x0, options);
u = input_2d(0,0,'allout');

%% Make results at a sampled time and add noise

Ts = 0.1;
sigma = 0.01;   % noise standard dev
mean =  0;  % noise mean

% get the state at sampled times
tq = (0 : Ts : tspan(2))';
yq = interp1(t, y, tq);
uq = interp1(t, u, tq);

% specify which states are "measured"
yobs = yq;

% inject noise to simulate measurement noise
noise = sigma .* randn(size(yobs)) + mean;
xq = yobs + noise;

% make

% put data into a struct
data = struct;
data.t = tq;
data.x = xq;
data.u = uq;



%% Control Input

function output = input_2d(t, x, flag)
    
    persistent u;
    
    unew = 5*(1 - exp(-0.1*t)) + sin(t);
    status = 0;
    
    % initialize
    if isempty(u)
        u = unew;
    end
    
    
    % determine output based on flag
    if isempty(flag)
        % Successful integration step! Generate a new value:
        u = [u; unew];
%         u_out = unew;   % no one will see this output
        output = status;
    elseif strcmp(flag, 'allout')
        output = u; 
    elseif strcmp(flag, 'clear')
        clear u;
        output = status;
    elseif strcmp(flag, 'pass2ode')
        output = u(end);
%         status = u_out;  % makes sure the ode reads in the right value
    end
    
end