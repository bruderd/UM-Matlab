function [ p, psi, N, xsym ] = def_polyLift( n, max_degree )
%def_polyLift: Defines the lifting function that lifts state variable x to
% space spanned by monomials with total degree less than or equal to
% max_degree.
%   e.g. 1 x1 x2 x1^2 x1x2 x2^2 ...

x = sym('x', [n, 1]);   % state variable x

% Number of mononials, i.e. dimenstion of p(x)
N = factorial(n + max_degree) / ( factorial(n) * factorial(max_degree) );

% matrix of exponents (N x n). Each row gives exponents for 1 monomial
exponents = zeros(1,n);
for i = 1:max_degree
   exponents = [exponents; partitions(i, ones(1,n))]; 
end

% create vector of orderd monomials (column vector)
for i = 1:N
    p(i,1) = get_monomial(x, exponents(i,:));
end

% define matrix of exponents: columns=monomial term, rows=dimension of x
psi = exponents';

% create the lifting function: x -> p(x)
matlabFunction(p, 'File', 'polyLift', 'Vars', {x});

% output the symbolic x variable
xsym = x;

end

