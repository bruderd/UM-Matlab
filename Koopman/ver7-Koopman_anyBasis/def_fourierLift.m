function params = def_fourierLift( params )
%def_fourierLift: Defines the lifting function that lifts state variable x to
% a fourier basis

[n, p, naug, maxDegree] = deal(params.n, params.p, params.naug, params.maxDegree);

x = sym('x', [n, 1]);   % state variable x
u = sym('u', [p, 1]);   % input vector u
xaug = [x; u];   % augmented state variable;

% Number of basis elements, i.e. dimenstion of p(x)
% N = 1 + 2*factorial(naug + maxDegree) / ( factorial(naug) * factorial(maxDegree) );
N = (1 + 2*maxDegree)^naug; 

% create sins of cosines of all the states
poop = sym( zeros(1+2*maxDegree , naug) );
for i = 1 : naug
    poop(1,i) = 1;
    for j = 1 : maxDegree
        poop(2*j,i)   = cos(2*pi*j*xaug(i));
        poop(2*j+1,i) = sin(2*pi*j*xaug(i)); 
    end
end

% define fourier basis vector
fourierBasis = poop(:,1);
for i = 2 : naug 
    fourierBasis = kron(fourierBasis, poop(:,i));
end

% create the lifting function: x -> p(x)
matlabFunction(fourierBasis, 'File', 'stateLift', 'Vars', {x, u});

%% define derivative of lifted state with respect to x

dlift = jacobian(fourierBasis,x);
matlabFunction(dlift, 'File', 'jacobianLift', 'Vars', {x,u});

%% output variables  
params.Basis = fourierBasis;    % symbolic vector of basis functions, p(x)
params.jacobianBasis = dlift;   % symbolic jacobian of the vector of basis monomials
params.N = N;   % dimension of polyBasis
params.x = x;   % symbolic state variable
params.u = u;   % symbolic input variable
params.xaug = xaug; % symbolic augmented state variable

end