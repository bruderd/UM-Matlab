function error = calcError( x, xdes, params )
%calcError: Calculate the error
%   Detailed explanation goes here

error = x - xdes;

%% More complicated version for full PID control
% error = struct;
% % saves all the old error values, and appends extra column
% error.p = [error_old.p, zeros(6,1)];
% error.i = [error_old.i, zeros(6,1)];
% error.k = [error_old.d, zeros(6,1)];
% 
% % caclulate the error
% error.p(:,k) = x - xdes; % current error, error_k
% error.i(:,k) = error_old.i(:,k-1) + error.p(:,k); % sum of current error and all previous errors since desired point was set, sum_0^k{ error_i }
% error.d(:,k) = error.p(:,k) - error_old.p(:,k-1); % difference between last two error measurements, (error_k - error_(k-1))

end

