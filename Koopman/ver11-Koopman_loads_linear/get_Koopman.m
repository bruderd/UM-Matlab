function [ U , koopData ] = get_Koopman( snapshotPairs , params )
%get_KoopmanConstGen: Find the best possible koopman operator given
%snapshot pairs using constraint generation to deal with large data sets.
%   Detailed explanation goes here

disp('Finding Koopman operator approximation...');

%% Extract snapshot pairs

% [x,y,u] = deal(snapshotPairs.x, snapshotPairs.y, snapshotPairs.u);
[x,y,u,w] = deal(snapshotPairs.zeta_x, snapshotPairs.zeta_y, snapshotPairs.u, snapshotPairs.w);

%% Build matrices

[n, p, nw] = deal(params.n, params.p, params.nw);

Px = zeros(length(x), params.Np + params.nw);
Py = zeros(length(x), params.Np + params.nw);
for i = 1:length(x)
    psix = stateLift( x(i,:)' )';
    psiy = stateLift( y(i,:)' )';
    Px(i,:) = [ psix , u(i,:) , w(i,:) ];
    Py(i,:) = [ psiy , u(i,:) , w(i,:) ];     % exclude u from Py (could also use same u as Px
end

K = size(Px,1);
Np = size(Px,2);

%% Store useful data that can be used outside this function
koopData.Px = Px( : , 1 : params.N );   % only want state, not input or load
koopData.Py = Py( : , 1 : params.N );
koopData.x = snapshotPairs.x;
koopData.u = u;
koopData.w = w;
koopData.zeta_x = snapshotPairs.zeta_x;

%% Solve for inital Koopman Operator wish subset of data points

% Add the dimension of load to Np before passing into Uvec (local change)
params.Np = params.Np + params.nw;

% Call function that solves QP problem
Uvec = solve_KoopmanQP(Px, Py, params);

U = reshape(Uvec, [Np,Np]);

%% check how well U works at each point (optional)
% dif = Px * U - Py;
% dif_x = dif( : , 1 : params.n);

disp('Done');

end

