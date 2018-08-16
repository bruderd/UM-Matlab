function [u, freq] = get_validationInput( t, params )
%set_vlidationParams: Set the value of parameters need to run sysid experiments
%   This is a slightly randomized version of get_sysidInput


Tamp_min = params.amp.minPeriod;    % minimum period of sinusoid that changes amplitude u_i
famp_max = 2*pi*(1/Tamp_min);  % maximum frequency of sinusoid that changes amplitude of u_i
T_min = params.freq.minPeriod;      % minimum period at which to drive any actuator
f_max = 2*pi*(1/T_min);   % maximum frequency at which to drive any actuator 

% initializations
amp = zeros(params.p,1);
freq = zeros(params.p,1);
u = zeros(params.p,1);

% offsets between inputs
ampoffset = params.val.ampoffset;
amp2freq = params.val.amp2freq;
foffset = params.val.foffset;

% add a delay at the beginning to give time for mocap to initialize etc.
t = t - params.tstart;       

if t < 0
    u = zeros(params.p,1);
else
    for i = 1 : params.p
        fmulti = 1 - foffset*(i-1);    % frequency multiplier
        ampmulti = (2.5 / params.p) * ampoffset * (i);     % amplitude multiplier
        
        % set amplitude of each input
        amp(i) = -ampmulti * cos(fmulti * famp_max * t) + ampmulti;
        
        % set frequency of each input
        foo = fmulti * amp2freq * famp_max;
        bar = f_max/2;
        freq(i) = -bar * cos(foo*t) + bar;
        
        % set outputs
        funt = bar * ( t - sin(foo*t) / foo );
        u(i,1) = amp(i) * sin(funt) + amp(i);
    end
end

end

