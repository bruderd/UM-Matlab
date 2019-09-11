function out = koop2cts( U , snapshotPairs , params )
%koop2cts: Converts the Koopman operator representation of a dynamical 
% system into a continuous nonlinear state space model.
% 
% See "Nonlinear system identification of soft robot dynamics" by Bruder
% et. al.

%% Simulate and find Koopman operator from "measurements" (DNE)

[x,y] = deal(snapshotPairs.zeta_x, snapshotPairs.zeta_y);
% [x,y, L_scale, R_scale] = scale_snapshotPairs( snapshotPairs , params );    % scale the snapshot pairs to help with model fitting

%% Calculate the infiniesimal generator as funtion of coeffients, and from data (DNE)
disp('Calculating infinitesimal generator...')
Ldata = get_Ldata(U, params);   % infinitesimal generator from data
disp('Done.')

%% solve for the coefficients, i.e. Eq. (18) from Mauroy and Gonclaves (DNE)

% create function that calculates the jacobian of the lifted state wrt. zeta
% (specifically only works for one delay)
params.zeta = [ params.x ; params.xd ; params.ud ];
params.Basis_aug = [ params.Basis ; params.u ];
jBasis = jacobian( params.Basis_aug , params.zeta );
params.jacobianLift = matlabFunction( jBasis , 'Vars' , { params.zeta } );

% matrix of coefficents of monomials
w = calc_W( Ldata , x , params );

% w = L_scale * w * R_scale;    % scale the coefficients back up so that they can explain dynamics of real model 

% dynamics (gives symbolic expression in terms of state and input)
vf2 = w * params.Basis_aug;

% create a function that evaluates continuous nonlinear dynamics
vf_cts = matlabFunction(vf2, 'File', 'vf_koopman', 'Vars', {params.zeta, params.u});


%% Define outputs
out.U       = U;            % koopman operator
out.Ldata   = Ldata;        % inf. generator from data
out.w       = w;            % matrix of coefficients of polybasis
out.vf      = vf2;          % symbolic vector field
out.vf_fh      = vf_cts;       % function handle for dynamics of sysid'd sys.

end
