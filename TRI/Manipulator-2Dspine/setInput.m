function u = setInput(t, params)
%Sets the input u(t) = P(t)
%   Detailed explanation goes here

u = zeros(params.n,1);
u(2) = (1-exp(-t))*1e12;

end

