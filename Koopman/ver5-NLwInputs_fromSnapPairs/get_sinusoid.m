function [ sinusoid ] = get_sinusoid( x, multipliers, params )
%get_sinusoid: builds a sinusoid from symbolic vector x and a vector of
%frequency multipliers
%   e.g. x = [x1 x2]; exponents = [1 2]; =>  sinusoid = sin(x1) * sin(2*x2)

n = length(x);
interval = params.interval;

sinusoid = cos( 2*pi * multipliers(1) * ( (x(1) - interval(1)) / (interval(2) - interval(1)) ) );
for i = 2:n
    sinusoid = sinusoid * cos( 2*pi * multipliers(i) * ( (x(i) - interval(1)) / (interval(2) - interval(1)) ) );
end

end