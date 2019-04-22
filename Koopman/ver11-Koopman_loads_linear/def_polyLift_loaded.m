function params = def_polyLift_loaded( params )
%def_polyLift_loaded: Defines the lifting function that lifts state variable x to
% space spanned by monomials with total degree less than or equal to
% max_degree.
%   e.g. 1 x1 x2 x1^2 x1x2 x2^2 ...

[n, p, nzeta, maxDegree] = deal(params.n, params.p, params.nzeta, params.maxDegree);

x = sym('x', [n, 1]);   % state variable x
xd = sym('xd', [params.nd * params.n, 1]);   % state delays i.e. for 2 delays: [x_i-1, x_i-2]'
ud = sym('ud', [params.nd * params.p, 1]);   % input delays i.e. for 2 delays: [u_i-1, u_i-2]'
w = sym('w', [params.nw , 1]);  % the load(s)
zeta = [x ; xd; ud];    % state variable with delays
u = sym('u', [p, 1]);   % input vector u

% matrix of exponents (N x naug). Each row gives exponents for 1 monomial
exponents = [];
for i = 1:maxDegree
   exponents = [exponents; partitions(i, ones(1,nzeta))]; 
end
exponents = [exponents ; zeros(1,nzeta)];   % put constant at end of basis so state can be the first nzeta elements

% create vector of orderd monomials (column vector)
Nnl = factorial(nzeta + maxDegree) / ( factorial(nzeta) * factorial(maxDegree) );
for i = 1:Nnl
    polyBasis(i,1) = get_monomial(zeta, exponents(i,:));
end

% create unloaded basis matrix
Wbasis = kron( eye(params.nw+1) , polyBasis );

% multiply the monomials by the load and append it to the basis vector
polyBasis = Wbasis * [ 1 ; w ];

% Number of mononials, i.e. dimenstion of the basis
N = size( polyBasis, 1);

% define matrix of exponents: columns=monomial term, rows=dimension of x
psi = exponents';

% create the lifting function: x -> p(x)
matlabFunction(polyBasis, 'File', 'stateLift', 'Vars', {zeta, w});

% create  a lifting function for the unloaded basis matrix
matlabFunction(Wbasis , 'File', 'Wlift' , 'Vars', {zeta});

%% define derivative of lifted state with respect to x

dlift = jacobian(polyBasis,x);
matlabFunction(dlift, 'File', 'jacobianLift', 'Vars', {zeta, w});

%% output variables  
params.Basis = polyBasis;    % symbolic vector of basis monomials, p(x)
params.jacobianBasis = dlift;
params.N = N;   % dimension of polyBasis (including the state itself)
params.Np = N + p; % dimension of the lifted state
params.psi = psi;   % monomial exponent index function
params.x = x;   % symbolic state variable
params.u = u;   % symbolic input variable
params.xd = xd; % symbolic state delays
params.ud = ud; % symbolic input delays
params.w = w;   % symbolic load variable
params.Wbasis = Wbasis; % symbolic unloaded basis matrix

end