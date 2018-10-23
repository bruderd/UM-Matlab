function params = def_fourierLift( params )
%def_fourierLift: Defines the lifting function that lifts state variable x to
% a fourier basis

[n, p, naug, maxDegree] = deal(params.n, params.p, params.naug, params.maxDegree);

x = sym('x', [n, 1]);   % state variable x
u = sym('u', [p, 1]);   % input vector u
% xaug = [x; u];   % augmented state variable;

% Number of basis elements, i.e. dimenstion of p(x)
% N = 1 + 2*factorial(naug + maxDegree) / ( factorial(naug) * factorial(maxDegree) );
N = n + (1 + 2*maxDegree)^n; 

% create sins of cosines of all the states
poop = sym( zeros(1+2*maxDegree , n) );
for i = 1 : n
    poop(1,i) = 1;
    for j = 1 : maxDegree
        poop(2*j,i)   = cos(2*pi*j*x(i));
        poop(2*j+1,i) = sin(2*pi*j*x(i)); 
    end
end

% define fourier basis vector
fourierBasis = poop(:,1);
for i = 2 : n 
    fourierBasis = kron(fourierBasis, poop(:,i));
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
params.N = N;   % dimension of basis (including the state itself)
params.Np = N + p;  % dimension of the lifted state
params.x = x;   % symbolic state variable
params.u = u;   % symbolic input variable
% params.xaug = xaug; % symbolic augmented state variable

end