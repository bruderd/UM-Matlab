function u = setInput( t, params )
%setInput: Specify the pressure input, u(t), for the system here
%   Detailed explanation goes here

p = params.p;
n = params.n;

u = ones(sum(n),1);

% u(2) = 1;
% u(3) = 0;

end

