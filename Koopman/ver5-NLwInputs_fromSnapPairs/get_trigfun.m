function [ sinusoid ] = get_trigfun( x, exponents, multipliers, params )
%get_sinusoid: builds a sinusoid from symbolic vector x and a vector of
%frequency multipliers
%   e.g. x = [x1 x2]; exponents = [1 2]; =>  sinusoid = sin(x1) * sin(2*x2)

n = length(x);
interval = params.interval;

% xshift = 2 .* pi .* (x - interval(1)) ./ (interval(2) - interval(1));  % shift each element of x into interval;
xshift = x;

sinusoid = 1;   % initialization
for i = 1 : n
   j = 2*i - 1;
   sinusoid = sinusoid * cos(multipliers(j) * xshift(i))^exponents(j) * sin(multipliers(j+1) * xshift(i))^exponents(j+1);
end

end