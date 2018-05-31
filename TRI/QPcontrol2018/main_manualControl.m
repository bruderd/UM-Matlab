% main_fbcontrol

error = struct; % initialize the error as a struct

%% Accept user input for desired points

%% Get state from Jihong's code and calculate the error, its integral and its derivative

error.p = ; % current error, error_k
error.i = ; % sum of current error and all previous errors since desired point was set, sum_0^k{ error_i }
error.d = ; % difference between last two error measurements, (error_k - error_(k-1))

%% Solve the QP for the desired pressure

%% Send commanded pressure values to Labview
%   can convert to voltage signals within the labview script

%% Wait specified amound of time

%% repeat until termination condition met