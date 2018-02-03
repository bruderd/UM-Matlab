function u = setInput(t, params)
%Sets the input u(t) = P(t)
%   Detailed explanation goes here

u = zeros(params.n,1);
u(1) = (1-exp(-t))*1e6;

end

