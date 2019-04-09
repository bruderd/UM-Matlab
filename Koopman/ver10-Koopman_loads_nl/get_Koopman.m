function [ U , koopData ] = get_Koopman( snapshotPairs , params )
%get_KoopmanConstGen: Find the best possible koopman operator given
%snapshot pairs using constraint generation to deal with large data sets.
%   Detailed explanation goes here

disp('Finding Koopman operator approximation...');

%% Extract snapshot pairs

% [x,y,u] = deal(snapshotPairs.x, snapshotPairs.y, snapshotPairs.u);
[x,y,u] = deal(snapshotPairs.zeta_x, snapshotPairs.zeta_y, snapshotPairs.u);

%% Build matrices

[n, p] = deal(params.n, params.p);

Px = zeros(length(x), params.Np);
Py = zeros(length(x), params.Np);
for i = 1:length(x)
    psix = stateLift( x(i,:)' )';
    psiy = stateLift( y(i,:)' )';
    Px(i,:) = [ psix , u(i,:) ];
    Py(i,:) = [ psiy , zeros(1,p) ];     % exclude u from Py (could also use same u as Px
end

K = size(Px,1);
Np = size(Px,2);

%% Store useful data that can be used outside this function
koopData.Px = Px( : , 1 : params.N );   % only want state, not input
koopData.Py = Py( : , 1 : params.N );
koopData.x = snapshotPairs.x;
koopData.u = u;
koopData.zeta_x = snapshotPairs.zeta_x;

%% Solve for inital Koopman Operator wish subset of data points

% Call function that solves QP problem
Uvec = solve_KoopmanQP(Px, Py, params);

U = reshape(Uvec, [Np,Np]);

%% check how well U works at each point (optional)
% dif = Px * U - Py;
% dif_x = dif( : , 1 : params.n);

disp('Done');

end

