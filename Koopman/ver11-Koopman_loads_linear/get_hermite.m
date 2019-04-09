function [ hermite ] = get_hermite( x , exponents )
%get_hermite: builds a monomial (product of hermites) from symbolic vector x and a vector of
%exponents
%   e.g. x = [x1 x2]; exponents = [1 2]; =>  monomial = x1^1 * x2^2

n = length(x);
m = size(exponents,2);  % width of the exponents matrix

% hermite = x(1)^exponents(1);
% for i = 2:n
%     hermite = hermite * x(i)^exponents(i);
% end

hermite = 1;
for i = 1:n
        hermite = hermite * hermiteH( exponents(i) , x(i) );
end

end