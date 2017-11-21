function u = setInput(t, params)
%Sets the input u(t) = P(t)
%   Detailed explanation goes here

u = zeros(params.n,1);
u(1) = exp(t-10)*100;

end

