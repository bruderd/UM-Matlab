function mpc = set_mpc( model , mpcParams )
%set_mpc: Creates the matrices for MPC optimal control problem
%   Detailed explanation goes here

%% define cost function matrices
% Cost function is defined: U'HU + ( z0'G + Yr'D )U

% A
N = size(model.A,1);
A = sparse( N*(mpcParams.Np+1) , N );
for i = 0 : mpcParams.Np
    A( (N*i + 1) : N*(i+1) , : ) = model.A^i ;
end

% B
Bheight = N*(mpcParams.Np+1);
Bcolwidth = size(model.B,2);
Bcol = sparse( Bheight , Bcolwidth );    % first column of B matrix
for i = 1 : mpcParams.Np
    Bcol( (N*i + 1) : N*(i+1) , : ) = model.A^(i-1) * model.B ;
end

Lshift = spdiags( ones( N*mpcParams.Np , 1 ) , -N , N*(mpcParams.Np+1) , N*(mpcParams.Np+1) );    % lower shift operator

Bwidth = size(model.B,2)*(mpcParams.Np);    % total columns in B matrix
Bblkwidth = mpcParams.Np;   % how many Bcol blocks wide B matrix is
B = spalloc( Bheight , Bwidth , floor(Bheight * Bwidth / 2) ); % initialze sparse B matrix
B(: , 1:Bcolwidth) = Bcol;
for i = 2 : Bblkwidth
    B(: , (i-1)*Bcolwidth+1 : i*Bcolwidth) = Lshift * B(: , (i-2)*Bcolwidth+1 : (i-1)*Bcolwidth);
end

% C
C = kron( speye(mpcParams.Np+1) , model.C);

% Q
Q = kron( speye(mpcParams.Np+1) , eye(model.params.ny) * 1); % error magnitude penalty

% R
R = kron( speye(mpcParams.Np) , eye(model.params.p) * 0.0001);  % input magnitude penalty

% H, G, D
H = B' * C' * Q * C * B + R;
G = 2 * A' * C' * Q * C * B;
D = -2 * Q * C * B;

%% define constraint matrices 
% Constraint equation is defined: LU + Mz0 <= c

nc = mpcParams.nc;  % number of constraints (for simply bounded inputs, nc = 2*model.params.p)
Np = mpcParams.Np;     % steps in horizon
umin = mpcParams.umin * min( model.params.uScaleFactor );
umax = mpcParams.umax * max( model.params.uScaleFactor );

% F : for input constraints
Fi = [ -speye(model.params.p) ; speye(model.params.p) ];    % diagonal element of F, for bounded inputs
F = sparse( nc * (mpcParams.Np+1) , size(B,2) );  % no constraints, all zeros
F( 1:nc*Np , 1:Np*model.params.p ) = kron( speye(Np) , Fi );     % fill in nonzeros

% E : for state constraints
E = sparse( nc * (mpcParams.Np+1) , size(B,1) );  % no constraints, all zeros

% c : the right side of constraint equation
ci = [ ones(model.params.p,1) * umin ; ones(model.params.p,1) * umax ]; % [umin, umax]' sized appropriately (since umin/max are scalar)
c = sparse( nc * (Np+1) , 1);    % sparse initialization of right side of constraint equation
c(1 : nc*Np) = kron( ones( Np , 1 ) , ci );     % fill in nonzeros

% L , M
L = F + E*B;
M = E*A;

%% save matrices in struct
mpc.H = H; mpc.G = G; mpc.D = D; mpc.L = L; mpc.M = M; mpc.c = c;

mpc.params = mpcParams; % also include the mpc parameters

end

