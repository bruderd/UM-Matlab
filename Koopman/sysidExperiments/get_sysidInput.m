function u = get_sysidInput( t, params )
%get_sysidInput: calculate the input signal (in volts) to the TR pressure
%regulators at each timestep. Meant to sweep over most of the control space
%   Detailed explanation goes here

Tamp_min = params.amp.minPeriod;    % minimum period of sinusoid that changes amplitude u_i
famp_max = 2*pi*(1/Tamp_min);  % maximum frequency of sinusoid that changes amplitude of u_i
T_min = params.freq.minPeriod;      % minimum period at which to drive any actuator
f_max = 2*pi*(1/T_min);   % maximum frequency at which to drive any actuator 

% initializations
amp = zeros(params.p,1);
freq = zeros(params.p,1);
u = zeros(params.p,1);

% add a delay at the beginning to give time for mocap to initialize etc.
t = t - params.tstart;       

if t < 0
    u = zeros(params.p,1);
else
    for i = 1 : params.p
        multi = 1 - (1/(2.5*params.p))*(i-1);    % frequency multiplier
        
        % set amplitude of each input
        amp(i) = -2.5 * cos(multi * famp_max * t) + 2.5;
        
        % set frequency of each input
        freq(i) = -f_max * cos(multi * 5 * famp_max*t) + f_max;
        
        % set outputs
        u(i,1) = amp(i) * sin(freq(i) * t) + amp(i);
    end
end


end

