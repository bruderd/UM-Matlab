function main_manualControl( filename, params )
% main_fbcontrol: control the end effector postion with user inputs as the desired states. 

error = struct; % initialize the error as a struct

tic;    % begin the timer so that data can be recorded as a timeseries
k = 0;  % initialize step counter

while(true) % <--- GIVE THIS A TERMINATION CONDITION
% Accept user input for desired points




% Get state from Jihong's code and calculate the error, its integral and its derivative

xmeas = getState(); % read in the state of the end effector from Jihong's C++ code using TCP/IP socket
tk = toc;   % current time
csvwrite(filename, [tk, xmeas'], k, 0); % store the measured state at each time step in a csvfile

error.p = ; % current error, error_k
error.i = ; % sum of current error and all previous errors since desired point was set, sum_0^k{ error_i }
error.d = ; % difference between last two error measurements, (error_k - error_(k-1))





% Solve the QP for the desired pressure





% Send commanded pressure values to Labview
%   can convert to voltage signals within the labview script





% Wait specified amound of time

k = k + 1;  % increase step counter
end

%% repeat until termination condition met