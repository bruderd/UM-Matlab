function params = set_sysidParams( filename )
%set_sysidParams: Set the value of parameters need to run sysid experiments
%   Detailed explanation goes here

params = struct;
amp = struct;   % for parameters related to amplitude of inputs
freq = struct;  % for parameters related to frequency of inputs
val = struct;   % for parameters rlated to validation

% Amplitude Sinusoid related parameters
amp.minPeriod = 240;

% Frequency Sinusoid related parameters
freq.minPeriod = 2;

% Validation parameters
val.amp2freq = 5 * rand;    % relates freq of amp sinusoid to signal sinusoid
val.foffset = rand;     % relates the frequencies between inputs
val.ampoffset = rand;   % relates the max amplitude of each signal

% Other parameters
params.p = 3;       % number of inputs
params.tstart = 10;     % wait this long before nonzero control inputs
params.amp = amp;
params.freq = freq;
params.val = val;


% Save parameters to file
filepath = ['paramFiles', filesep, filename, '.mat'];
save(filepath, 'params');



end

