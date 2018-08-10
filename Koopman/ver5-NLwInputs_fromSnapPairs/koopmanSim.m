function [ t, x ] = koopmanSim( U, x0, params )
%koopmanSim: Flow the system forward using the Koopman operator
%   Detailed explanation goes here

t = (0 : params.Ts : params.duration)';

% selection matrix to pull states out of lifted state
select = [zeros(params.n,1), [eye(params.n/2) ; zeros(params.n/2)], zeros(params.n, params.N - params.n/2 - 1)];

z = zeros(params.N, length(t));
x = zeros(length(t), params.n);
z(:,1) = polyLift(x0, 0);   % assume input is zero
x(1,:) = x0';
for i = 2 : length(t)
   z(:,i) = U' * z(:,i-1);
   xt = select * z(:,i);
   x(i, :) = xt';
end

end

