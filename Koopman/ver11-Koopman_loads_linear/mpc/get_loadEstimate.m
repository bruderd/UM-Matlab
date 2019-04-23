function [ load , Whor ] = get_loadEstimate( ynow , Wold, mpcData , mpc , model )
%get_loadEstimate: Estimates the load based on model and measurements
%   Detailed explanation goes here

Nw = mpc.params.Nw;
nd = model.params.nd;   % number of delays in model
scale_y = diag( model.params.xScaleFactor );
scale_u = diag( model.params.uScaleFactor );

% scale the measurements
ynow_sc = diag(model.params.xScaleFactor) * ynow;
if size(mpcData.Y,1) > nd + 1
    yd_sc = mpcData.Y( end-1-nd+1 : end-1 , :) * scale_y; % only because ynow has already been saved in mpcData
    ud_sc = mpcData.U( end-nd+1 : end , : ) * scale_u;
else
    yd_sc = zeros( nd , model.params.ny );
    ud_sc = zeros( nd , model.params.p );
end
ylast = reshape( fliplr( yd_sc' ) , [nd*model.params.ny,1] ); 
ulast = reshape( fliplr( ud_sc' ) , [nd*model.params.p,1] );

% construct zeta
if nd == 0
    zeta = ynow_sc;
else
    zeta = [ ynow_sc ;...
             ylast ;...
             ulast ];
end
    
%% lift the state (must first identify lifting function corresponding to current model)
cd( ['..' , filesep , 'liftingFunctions'] );
Wlift = str2func( [ 'Wlift_' , model.params.systemName ] );
cd(['..' , filesep , 'mpc']);
Wnow = Wlift(zeta);   % evaluate the state dependent load matrix

%% get estimate of the load using least squares over past horizon (on xhat)
Whor = [ Wold(model.params.N+1 : end , :) ; Wnow ]; % update the concatentaion of lifted load matrices

lendata = size(mpcData.Y , 1);  % number of recorded steps
if lendata > Nw
    Yhor = mpcData.Y( end - Nw + 1 : end , : ) * scale_y;
    Uhor = mpcData.U( end - Nw + 1 : end , : ) * scale_u;
else
    Yhor = [ zeros( Nw - lendata , model.params.ny ) ; mpcData.Y( end - lendata + 1 : end , : ) ] * scale_y;
    Uhor = [ zeros( Nw - lendata , model.params.p ) ; mpcData.U( end - lendata + 1 : end , : ) ] * scale_u;  
end
Yhor_vec = reshape( fliplr( Yhor' ) , [Nw*model.params.ny , 1] ); 
Uhor_vec = reshape( fliplr( Uhor' ) , [Nw*model.params.p , 1] );

Clsqlin = mpc.west.CAstack * Whor;
dlsqlin = Yhor_vec - mpc.west.CBstack * Uhor_vec;
lb = zeros(model.params.nw+1,1);  % load cannot be negative
ub = [ Inf ; model.params.scale * ones(model.params.nw,1) ]; % load should be bounded by maximum observed
sol = lsqlin( Clsqlin , dlsqlin , [] , [] , mpc.west.Aeq , mpc.west.beq , lb , ub );  % solve for what using constrained least squares solver
load = sol(2:end);  % the firt element must always be one so we skip over it


end

