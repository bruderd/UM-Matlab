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

% Q: Error magnitude penalty (CAN MODIFY THINGS HERE)
Q = kron( speye(mpcParams.Np+1) , eye(model.params.ny) * 0.1); % error magnitude penalty (running cost) (default 0.1)
% Q = sparse( (mpcParams.Np+1)*model.params.ny , (mpcParams.Np+1)*model.params.ny );  % all zeros
Q(end-model.params.ny+1 : end , end-model.params.ny+1 : end) = eye(model.params.ny) * 100;    % (terminal cost) (default 100)

% R: Input magnitude penalty
R = kron( speye(mpcParams.Np) , eye(model.params.p) * 0.5e-2 );  % input magnitude penalty (for flaccy use 0.5e-2)

% H, G, D
H = B' * C' * Q * C * B + R;
G = 2 * A' * C' * Q * C * B;
D = -2 * Q * C * B;

%% define constraint matrices 
% Constraint equation is defined: LU + Mz0 <= c

nc = mpcParams.nc;  % number of constraints (for simply bounded inputs, nc = 2*model.params.p)
Np = mpcParams.Np;     % steps in horizon
umin = mpcParams.umin * max( model.params.uScaleFactor );
umax = mpcParams.umax * min( model.params.uScaleFactor );

%% F : for input constraints

Fi = [ -speye(model.params.p) ; speye(model.params.p) ];    % diagonal element of F, for bounded inputs
F = sparse( nc * (mpcParams.Np+1) , size(B,2) );  % no constraints, all zeros
F( 1:nc*Np , 1:Np*model.params.p ) = kron( speye(Np) , Fi );     % fill in nonzeros

% %Add on input slope constraint (see pages 181-182 in lab notebook)
Fslope_i = speye(model.params.p);
Fslope_neg = [ kron( speye(Np-1) , -Fslope_i ) , sparse( model.params.p * (Np-1) , model.params.p ) ];
Fslope_pos = [ sparse( model.params.p * (Np-1) , model.params.p ) , kron( speye(Np-1) , Fslope_i ) ];
Fslope_top = Fslope_neg + Fslope_pos;
Fslope = [ Fslope_top ; -Fslope_top];
F = [ F ; Fslope ];

% % %Add on '2 actuator only' constraint
% F2only = kron( speye(Np) , ones(1,model.params.p ) );
% F = [ F ; F2only ];

% % Add on input smoothness constraint (see page 6 in new lab notebook)
Fsmooth_i = speye(model.params.p);
Fsmooth_lI = [ kron( speye(Np-2) , Fsmooth_i ) , sparse( model.params.p * (Np-2) , 2 * model.params.p ) ];
Fsmooth_2I = [ sparse( model.params.p * (Np-2) , model.params.p ) , kron( speye(Np-2) , -2*Fslope_i ) , sparse( model.params.p * (Np-2) , model.params.p ) ];
Fsmooth_rI = [ sparse( model.params.p * (Np-2) , 2 * model.params.p ) , kron( speye(Np-2) , Fslope_i ) ];
Fsmooth_top = Fsmooth_lI + Fsmooth_2I + Fsmooth_rI;
Fsmooth = [ Fsmooth_top ; -Fsmooth_top ];
F = [ F ; Fsmooth ];

%% E : for state constraints

E = sparse( nc * (mpcParams.Np+1) , size(B,1) );  % no constraints, all zeros

% %Add dimensions to make size consistent with input slope constraints
E = [ E ; sparse( 2 * model.params.p * (mpcParams.Np-1) , size(B,1) ) ];

% % %Add dimensions to make size consistent with input 2 actuator only constraint
% E = [ E ; sparse( Np , size(E,2) ) ];

% % Add dimensions to make size consistent with input smoothness constraint
E = [ E ; sparse( 2 * model.params.p * (mpcParams.Np-2) , size(B,1) ) ];


%% c : the right side of constraint equation (HERE YOU CAN MODIFY SLOPE AND SMOOTHNESS)

ci = [ -ones(model.params.p,1) * umin ; ones(model.params.p,1) * umax ]; % [umin, umax]' sized appropriately (since umin/max are scalar)
c = sparse( nc * (Np+1) , 1);    % sparse initialization of right side of constraint equation
c(1 : nc*Np) = kron( ones( Np , 1 ) , ci );     % fill in nonzeros

% %Add RHS of input slope constraints (limit imputs to the 30 second ramp)
slope_lim = ( (umax - umin)/30 ) * model.params.Ts;
cslope_top = slope_lim * ones( model.params.p * (Np-1) , 1 );
cslope = [ cslope_top ; cslope_top ];
c = [ c ; cslope ];

% % %Add RHS of '2 actuator only constraint'
% c2only = ones(Np,1) * 2*umax;
% c = [ c ; c2only ];

% %Add RHS of smoothness constraint
smooth_lim = model.params.Ts^2 * 0.04;      % for real system use model.params.Ts^2 * 0.04
csmooth = smooth_lim * ones( size(Fsmooth,1) ,1);
c = [ c ; csmooth ];

%% L , M
L = F + E*B;
M = E*A;

%% Matrices for the load estimator

% build matrices for estimating the load
whor = mpcParams.Nw;    % number of time steps to be considered in load estimate equation
mpc.west.Astack = kron(eye(whor) , model.A);
mpc.west.Bstack = kron(eye(whor) , model.B);
mpc.west.CAstack = kron(eye(whor) , model.Cy * model.A);
mpc.west.CBstack = kron(eye(whor) , model.Cy * model.B);
mpc.west.Aeq = blkdiag( 1 , zeros(model.params.nw , model.params.nw) );
mpc.west.beq = [ 1 ; zeros(model.params.nw,1) ];

%% save matrices in struct
mpc.H = H; mpc.G = G; mpc.D = D; mpc.L = L; mpc.M = M; mpc.c = c;

% other matrices for debugging
mpc.A = A; mpc.B = B; mpc.C = C; mpc.Q = Q; mpc.R = R;

mpc.params = mpcParams; % also include the mpc parameters

end

