function params = def_trigLift( params )
%def_sinLift: Defines the lifting function that lifts state variable x to
% space spanned by sines/cosines with frequency multiple less than or equal
% to maxDegree.
%   e.g. 1, sin(x1), cos(x1), ...

[n, p, naug, maxDegree] = deal(params.n, params.p, params.naug, params.maxDegree);
interval = params.interval;

x = sym('x', [n, 1]);   % state variable x
u = sym('u', [p, 1]);   % input vector u
xaug = [x; u];   % augmented state variable;


% selection matrix. Each row gives index for which trig terms are included
exponents = permn([0,1], 2*naug);
Nexp = size(exponents, 1);

% multiplier matrix. Each row gives the frequency multiplier of each trig term
multipliers = permn(1:maxDegree, 2*naug);
Nmult = size(multipliers, 1);

% create vector of orderd sinusoids (column vector)
for i = 1 : Nexp
    for j = 1 : Nmult 
        index = (i-1) * Nmult + j;
        trigBasis(index,1) = get_trigfun(xaug, exponents(i,:), multipliers(j,:), params);
    end
end

% define matrix of exponents: columns=monomial term, rows=dimension of x
psi = exponents';

% create the lifting function: x -> p(x)
matlabFunction(trigBasis, 'File', 'trigLift', 'Vars', {x, u});

% Number of mononials, i.e. dimenstion of p(x)
N = size(trigBasis, 1);

% Number of monomials in basis set for observables to be mapped through Lkj
N1 = factorial(naug + params.m1) / ( factorial(naug) * factorial(params.m1) );

% output variables  
params.trigBasis = trigBasis;    % symbolic vector of basis monomials, p(x)
params.N = N;   % dimension of sinBasis
params.N1 = N1; % dimension of basis of observables to be mapped through Lkj
params.psi = psi;   % monomial exponent index function
params.x = x;   % symbolic state variable
params.u = u;   % symbolic input variable
params.xaug = xaug; % symbolic augmented state variable

end