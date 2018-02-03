function u = setInput( t, params )
%setInput: Specify the pressure input, u(t), for the system here
%   Detailed explanation goes here

p = params.p;
n = params.n;

u = zeros(sum(n),1);

u(1) = 0;
u(11) = 0; % * (1 - exp(-t));
u(3) = 0;
u(12) = 0; % * (1 - exp(-t));
% u(5) = -1 * (1 - exp(-t));
% u(7) = -1 * (1 - exp(-t));
u(1) = 10000;
u(3) = 10000;

end

