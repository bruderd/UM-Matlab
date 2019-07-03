function u = calc_u( t, params )
%Calculate the input to the system, which is a vector of pressures
%   Detailed explanation goes here

% u = [1 0 0 1]' * 1e6;

% u = [1 0 1 0]' * 1000*(sin(t) + 1);

% u = [1 0 1 0]' * 1000 * (sin(pi*t/5) + 1);

% u = [0.6 0 0.4 0.8]' * 1e8 * (sin(pi*t/5)*exp(-t*2e-1) + 1);

% u = [0.6 0 0.4 0.8]' * 1e8 * (1 - exp(-t * 2e0));

%% Do a sequence of steady control inputs

% mult = 3e6;
% if t < 2
%     u = [1 1 0 0]' * mult;
% elseif t < 4
%     u = [1 0 1 0]' * mult;
% elseif t < 6
%     u = [1 0 0 1]' * mult;
% elseif t < 8
%     u = [0 1 1 0]' * mult;
% elseif t < 10
%     u = [0 1 0 1]' * mult;
% elseif t <= 12
%     u = [0 0 1 1]' * mult;
% end

%% Drive each free with a sinusoid
% Pmax = 3e6;
% 
% cmd1 = (cos(t)*Pmax + Pmax)*0.5;
% cmd2 = (cos(t + pi/2)*Pmax + Pmax)*0.5;
% cmd3 = (cos(t + 2*pi/2)*Pmax + Pmax)*0.5;
% cmd4 = (cos(t + 3*pi/2)*Pmax + Pmax)*0.5;
% 
% u = [cmd1, cmd2, cmd3, cmd4]';


%% Control inputs to match labview (for TRI demo)

% Pmax = 2e6;
Pmax = params.Pmax;
flow = 10;
ti = mod(t,2);

% if t < 2
%     vi = sin(ti*pi/2);
%     u = [1 1 0 0]' * (Pmax - Pmax*exp(-vi*flow));
% elseif t < 4
%     vi = sin(ti*pi/2);
%     u = [1 0 1 0]' * (Pmax - Pmax*exp(-vi*flow));
% elseif t < 6
%     vi = sin(ti*pi/2);
%     u = [1 0 0 1]' * (Pmax - Pmax*exp(-vi*flow));
% elseif t < 8
%     vi = sin(ti*pi/2);
%     u = [0 1 1 0]' * (Pmax - Pmax*exp(-vi*flow));
% elseif t < 10
%     vi = sin(ti*pi/2);
%     u = [0 1 0 1]' * (Pmax - Pmax*exp(-vi*flow));
% elseif t <= 12
%     vi = sin(ti*pi/2);
%     u = [0 0 1 1]' * (Pmax - Pmax*exp(-vi*flow));
% end


% % Reordered to sync up with video of real system
% if t < 2
%     vi = sin(ti*pi/2);
%     u = [0 1 0 1]' * (Pmax - Pmax*exp(-vi*flow));
% elseif t < 4
%     vi = sin(ti*pi/2);
%     u = [1 0 0 1]' * (Pmax - Pmax*exp(-vi*flow));
% elseif t < 6
%     vi = sin(ti*pi/2);
%     u = [0 1 1 0]' * (Pmax - Pmax*exp(-vi*flow));
% elseif t < 8
%     vi = sin(ti*pi/2);
%     u = [1 0 1 0]' * (Pmax - Pmax*exp(-vi*flow));
% elseif t < 10
%     vi = sin(ti*pi/2);
%     u = [1 1 0 0]' * (Pmax - Pmax*exp(-vi*flow));
% elseif t <= 12
%     vi = sin(ti*pi/2);
%     u = [0 0 1 1]' * (Pmax - Pmax*exp(-vi*flow));
% end


% Reordered to sync up with video of real system (I mixed up tubes 2 and 3)
if t < 2
    vi = sin(ti*pi/2);
    u = [0 0 1 1]' * (Pmax - Pmax*exp(-vi*flow));
elseif t < 4
    vi = sin(ti*pi/2);
    u = [1 0 1 0]' * (Pmax - Pmax*exp(-vi*flow));
elseif t < 6
    vi = sin(ti*pi/2);
    u = [0 1 1 0]' * (Pmax - Pmax*exp(-vi*flow));
elseif t < 8
    vi = sin(ti*pi/2);
    u = [1 0 0 1]' * (Pmax - Pmax*exp(-vi*flow));
elseif t < 10
    vi = sin(ti*pi/2);
    u = [0 1 0 1]' * (Pmax - Pmax*exp(-vi*flow));
elseif t <= 12
    vi = sin(ti*pi/2);
    u = [1 1 0 0]' * (Pmax - Pmax*exp(-vi*flow));
end

end

