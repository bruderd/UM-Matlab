function out = koopmanSysid( snapshotPairs, params )
%koopmanSysid: Identifies system dynamics from snapshot pairs as described
% in Goncalves, Mauroy, "Linear identification of nonlinear systems: 
% A lifting technique based on the Koopman operator"
%   Detailed explanation goes here

%% Simulate and find Koopman operator from "measurements" (DNE)

[x,y] = deal(snapshotPairs.x, snapshotPairs.y);
% [x,y, L_scale, R_scale] = scale_snapshotPairs( snapshotPairs , params );    % scale the snapshot pairs to help with model fitting
U = get_Koopman(x,y, params);
% U = get_KoopmanRobust( x, y, params );


%% Calculate the infiniesimal generator as funtion of coeffients, and from data (DNE)
Ldata = get_Ldata(U, params);   % infinitesimal generator from data
Ldata_hat = Ldata(:, 1:params.N1);  % N x N1 version of Ldata (projection onto the polyBasis)
vecLdata = Ldata_hat(:);    % vectorized version of Ldata matrix

vecstackL = zeros(params.N*params.N1, params.N*params.n);
for k = 1:params.N
    for j = 1:params.n
        Lkj = get_Lkj(k,j,params);
        
        % convert all the Lkj's into vectors and stack them horizontally
        vecLkj = Lkj(:);
        vecstackL(:, (k-1)*params.n + j) = vecLkj;
    end
end

%% solve for the coefficients, i.e. Eq. (18) from Mauroy and Gonclaves (DNE)

W = pinv(vecstackL) * vecLdata;

% matrix of coefficents of monomials
w = reshape(W, [params.n, params.N]);
% w = L_scale * w * R_scale;    % scale the coefficients back up so that they can explain dynamics of real model

% dynamics (gives symbolic expression in terms of state and input)
vf2 = w * params.polyBasis; 
matlabFunction(vf2, 'File', 'vf_koopman', 'Vars', {params.x, params.u});


%% Define outputs
out.U       = U;            % koopman operator
out.Ldata   = Ldata;        % inf. generator from data
out.w       = w;            % matrix of coefficients of polybasis
out.vf      = @vf_koopman;  % function handle for dynamics of sysid'd sys.

end

