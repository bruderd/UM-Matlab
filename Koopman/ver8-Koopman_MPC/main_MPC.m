% main_MPC
%
%
%
%
%

%% load in system model
[model_file , model_path] = uigetfile('./models/');
model = load( [model_path , model_file] );

%% set MPC parameters
mpc = struct;

mpc.Ts      = model.params.Ts;      % sampling time must be the same as discrete model
mpc.tf      = 20;                     % total length of MPC simulation
mpc.tspan   = 0 : mpc.Ts : mpc.tf;  % time vector for MPC simulation
mpc.Np      = 100;                 % prediction horizon
mpc.x0      = zeros(model.params.n , 1);    % initial condition
mpc.nc      = 1;                % number of constraints


%% define reference trajectory

for i = 1 : length(mpc.tspan)
%     ref = [ sin( mpc.tspan(i) ) ; 0 ];
    ref = [ 1 ; 0 ];      % step input
    Yr( (i-1)*model.params.n+1 : i*model.params.n , 1 ) = ref;
end


%% define cost function matrices
% Cost function is defined: U'HU + ( z0'G + Yr'D )U

% A
N = size(model.A,1);
A = sparse( N*(mpc.Np+1) , N );
for i = 0 : mpc.Np
    A( (N*i + 1) : N*(i+1) , : ) = model.A^i ;
end

% B
Bheight = N*(mpc.Np+1);
Bcol = sparse( Bheight , size(model.B,2) );    % first column of B matrix
for i = 1 : mpc.Np
    Bcol( (N*i + 1) : N*(i+1) , : ) = model.A^(i-1) * model.B ;
end

Lshift = spdiags( ones( N*mpc.Np , 1 ) , -N , N*(mpc.Np+1) , N*(mpc.Np+1) );    % lower shift operator

Bwidth = size(model.B,2)*(mpc.Np);
B = spalloc( Bheight , Bwidth , floor(Bheight * Bwidth / 2) ); % initialze sparse B matrix
B(:,1) = Bcol;
for i = 2 : Bwidth
    B(:,i) = Lshift * B(:,i-1);
end

% C
C = kron( speye(mpc.Np+1) , model.C);

% Q
Q = kron( speye(mpc.Np+1) , eye(model.params.n)); % error magnitude penalty

% R
R = kron( speye(mpc.Np) , eye(model.params.p) * 0.01);  % input magnitude penalty

% H, G, D
H = B' * C' * Q * C * B + R;
G = 2 * A' * C' * Q * C * B;
D = -2 * Q * C * B;

%% define constraint matrices 

nc = mpc.nc;

% F
F = sparse( nc * (mpc.Np+1) , size(B,2) );  % no constraints, all zeros

% E
E = sparse( nc * (mpc.Np+1) , size(B,1) );  % no constraints, all zeros

% L , M
L = F + E*B;
M = E*A;

%% save matrices in struct
mpc.H = H; mpc.G = G; mpc.D = D; mpc.L = L; mpc.M = M;


%% run MPC simulation

% identify system dynamics function
[vf_file , vf_path] = uigetfile('./fakeSystems/simDynamics');
fileparsed = split( vf_file , '.' );
system_dynamics = str2func( fileparsed{1} );

addpath(vf_path);
[t, x] = ode45( @(t,x) system_dynamics(x, get_MPCinput(t,x,Yr,mpc,model)), mpc.tspan(1 : end-mpc.Np), mpc.x0 );
rmpath(vf_path);




% function to determine input at each time step
function unext = get_MPCinput( tnow , xnow , ref , mpc , model)

persistent k x u

know = floor( tnow / mpc.Ts );    % current timestep

% initialize values of persistent variables
if isempty(k)
    k = know;
    x = kron( ones(model.params.nd+1 , 1) , mpc.x0 );   % repeat initial condition for time steps that occur before k = 0
    u = kron( ones(model.params.nd+1 , 1) , zeros(model.params.p , 1) );    % assume input is zero before k = 0
end

% update values of persistent variables
if know == k(end)       % still at same timestep?
    unext = u(: , end);     % u takes same value 
else
    nd = model.params.nd;
    if nd == 0
        zeta = xnow;
    else
        zeta = [ xnow ; x( : , end-nd+1 : end ) ; u( : , end-nd+1 : end )];
    end
    z = stateLift(zeta);   % lift x
    Yr = ref( (know)*model.params.n+1 : ((know+1) + mpc.Np)*model.params.n , 1 ); % reference trajectory over prediction horizon
    H = 2 * mpc.H;
    f = ( z' * mpc.G + Yr' * mpc.D )';
    A = mpc.L;
    b = - mpc.M * z;
    U = quadprog( H , f , A , b );     % use H, G etc here...
    unext = U( 1 : model.params.p );    % the actual output of this function
    
    % store values of the state and input at the "beginning" of a timestep
    k = [k , know];     % vector of all previous timestamps
    x = [x , xnow];     % vector of all previous states
    u = [u , unext];
end



end




