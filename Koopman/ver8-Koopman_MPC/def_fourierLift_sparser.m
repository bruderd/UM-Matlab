function params = def_fourierLift_sparser( params )
%def_fourierLift: Defines the lifting function that lifts state variable x to
% a fourier basis

[n, p, naug, maxDegree] = deal(params.n, params.p, params.naug, params.maxDegree);

x = sym('x', [n, 1]);   % state variable x
u = sym('u', [p, 1]);   % input vector u
% xaug = [x; u];   % augmented state variable;

% matrix of exponents (N x naug). Each row gives exponents for 1 monomial
multipliers = zeros(1,2*n);
for i = 1:maxDegree
   multipliers = [multipliers; partitions(i, ones(1, 2*n))]; 
end

% Number of basis elements, i.e. dimenstion of p(x)
N = n + size(multipliers , 1);

% create vector of sines and cosines with multipliers
fourierBasis = sym('fourierBasis', [N-n,1]);
for i = 1:N-n
    fourierBasis(i,1) = get_sinusoid(x, multipliers(i,:));
end

% put the full state at the beginnig of the basis vector
fourierBasis = [x ; fourierBasis];

% create the lifting function: x -> p(x)
matlabFunction(fourierBasis, 'File', 'stateLift', 'Vars', {x});

%% define derivative of lifted state with respect to x

dlift = jacobian(fourierBasis,x);
matlabFunction(dlift, 'File', 'jacobianLift', 'Vars', {x});

%% output variables  
params.Basis = fourierBasis;    % symbolic vector of basis functions, p(x)
params.jacobianBasis = dlift;   % symbolic jacobian of the vector of basis monomials
params.N = N;   % dimension of basis
params.Np = N + p;  % dimension of lifted state
params.x = x;   % symbolic state variable
params.u = u;   % symbolic input variable
% params.xaug = xaug; % symbolic augmented state variable

end