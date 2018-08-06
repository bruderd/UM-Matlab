function params = def_sinLift( params )
%def_sinLift: Defines the lifting function that lifts state variable x to
% space spanned by sines with frequency multiple less than or equal
% to maxDegree.
%   e.g. 1, sin(x1), sin(x2), sin(2*x1), sin(x1)*sin(x2), sin(2*x2)  ...

[n, p, naug, maxDegree] = deal(params.n, params.p, params.naug, params.maxDegree);
interval = params.interval;

x = sym('x', [n, 1]);   % state variable x
u = sym('u', [p, 1]);   % input vector u
xaug = [x; u];   % augmented state variable;

% Number of mononials, i.e. dimenstion of p(x)
N = factorial(naug + maxDegree) / ( factorial(naug) * factorial(maxDegree) );

% Number of monomials in basis set for observables to be mapped through Lkj
N1 = factorial(naug + params.m1) / ( factorial(naug) * factorial(params.m1) );

% matrix of exponents (N x naug). Each row gives exponents for 1 monomial
exponents = zeros(1,naug);
for i = 1:maxDegree
   exponents = [exponents; partitions(i, ones(1,naug))]; 
end

% create vector of orderd sinusoids (column vector)
for i = 1:N
    sinBasis(i,1) = get_sinusoid(xaug, exponents(i,:), params);
end

% define matrix of exponents: columns=monomial term, rows=dimension of x
psi = exponents';

% create the lifting function: x -> p(x)
matlabFunction(sinBasis, 'File', 'sinLift', 'Vars', {x, u});

% output variables  
params.sinBasis = sinBasis;    % symbolic vector of basis monomials, p(x)
params.N = N;   % dimension of sinBasis
params.N1 = N1; % dimension of basis of observables to be mapped through Lkj
params.psi = psi;   % monomial exponent index function
params.x = x;   % symbolic state variable
params.u = u;   % symbolic input variable
params.xaug = xaug; % symbolic augmented state variable

end