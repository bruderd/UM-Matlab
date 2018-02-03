function x = setState( t, params )
%setInput: Specify the pressure input, u(t), for the system here
%   Detailed explanation goes here

p = params.p;
n = params.n;

x = zeros(3*p,1);

x(2) = 0.2761;

end