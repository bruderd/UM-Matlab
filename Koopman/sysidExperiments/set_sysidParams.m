function params = set_sysidParams( filename )
%set_sysidParams: Set the value of parameters need to run sysid experiments
%   Detailed explanation goes here

params = struct;
amp = struct;   % for parameters related to amplitude of inputs
freq = struct;  % for parameters related to frequency of inputs

% Amplitude Sinusoid related parameters
amp.minPeriod = 240;

% Frequency Sinusoid related parameters
freq.minPeriod = 2;

% Other parameters
params.p = 3;       % number of inputs
params.tstart = 10;     % wait this long before nonzero control inputs
params.amp = amp;
params.freq = freq;


% Save parameters to file
filepath = ['paramFiles', filesep, filename, '.mat'];
save(filepath, 'params');



end

