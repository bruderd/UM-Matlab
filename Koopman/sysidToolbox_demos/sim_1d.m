%% Simulate 1D system
clear;


tspan = [0 100];
x0 = 0;
options = odeset('OutputFcn',@(t,x,flag) input_1d(t, x, flag), 'refine', 1);

input_1d(0,0,'clear');
[t, y] = ode45(@(t,x) myode_1d(t, x, input_1d(t, x, 'pass2ode')), tspan, x0, options);
u = input_1d(0,0,'allout');


%% Control Input

function output = input_1d(t, x, flag)
    
    persistent u;
    
    unew = sin(t);
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