function params = def_polyLift( params )
%def_polyLift: Defines the lifting function that lifts state variable x to
% space spanned by monomials with total degree less than or equal to
% max_degree.
%   e.g. 1 x1 x2 x1^2 x1x2 x2^2 ...

[n, p, naug, maxDegree] = deal(params.n, params.p, params.naug, params.maxDegree);

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

% create vector of orderd monomials (column vector)
for i = 1:N
    polyBasis(i,1) = get_monomial(xaug, exponents(i,:));
end

% define matrix of exponents: columns=monomial term, rows=dimension of x
psi = exponents';

% create the lifting function: x -> p(x)
matlabFunction(polyBasis, 'File', 'polyLift', 'Vars', {x, u});

%% define derivative of lifted state with respect to x

dlift = jacobian(polyBasis,x);
matlabFunction(dlift, 'File', 'jacobianLift', 'Vars', {x,u});

%% output variables  
params.polyBasis = polyBasis;    % symbolic vector of basis monomials, p(x)
params.jacobianBasis = dlift;
params.N = N;   % dimension of polyBasis
params.N1 = N1; % dimension of basis of observables to be mapped through Lkj
params.psi = psi;   % monomial exponent index function
params.x = x;   % symbolic state variable
params.u = u;   % symbolic input variable
params.xaug = xaug; % symbolic augmented state variable

end