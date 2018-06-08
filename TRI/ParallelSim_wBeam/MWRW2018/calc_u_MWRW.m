function [u, error] = calc_u_MWRW( xin, t, params )
%Calculate the input to the system, which is a vector of pressures
%   Also calculates the error, which will be used for PID terms
%   Note: the input "xin" is the euler angles and their derivatives

% calculate full state
xeul = xin(1:3);
xcart = euler2cart(xeul, params);
x = [xcart; xeul];

%% Set the desired position

% Define the desired position
if t <= 12
    xeul_des = [0, 0, pi/8]';
    xcart_des = euler2cart(xeul_des, params);
    xdes = [xcart_des; xeul_des];
elseif t < 4
    xeul_des = [pi/8, 0, 0]';
    xcart_des = euler2cart(xeul_des, params);
    xdes = [xcart_des; xeul_des];
elseif t < 6
    xeul_des = [-pi/8, -pi/8, 0]';
    xcart_des = euler2cart(xeul_des, params);
    xdes = [xcart_des; xeul_des];
elseif t < 8
    xeul_des = [0, pi/8, pi/4]';
    xcart_des = euler2cart(xeul_des, params);
    xdes = [xcart_des; xeul_des];
elseif t < 10
    xeul_des = [0, pi/8, -pi/4]';
    xcart_des = euler2cart(xeul_des, params);
    xdes = [xcart_des; xeul_des];
elseif t <= 12
    xeul_des = [0, pi/4, 0]';
    xcart_des = euler2cart(xeul_des, params);
    xdes = [xcart_des; xeul_des];
end


%% Calculate the input pressure to reduce error

error = calcError(x, xdes, params);

fload = zeros(6,1); % no load for now, may add in weight of end effector later
[psol, exitflag] = calcPressure( x, error, params, fload);


% if infeasable set input equal to zero
if exitflag < 0
    u = zeros(4,1);
else
    u = psol(1:4,:) * 1.7; % put in an extra gain
end


end

