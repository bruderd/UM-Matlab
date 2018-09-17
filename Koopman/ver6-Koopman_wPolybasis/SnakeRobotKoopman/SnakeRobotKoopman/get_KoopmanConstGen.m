function U = get_KoopmanConstGen( x,y, params )
%get_KoopmanConstGen: Find the best possible koopman operator given
%snapshot pairs using constraint generation to deal with large data sets.
%   Detailed explanation goes here

[n, p] = deal(params.n, params.p);

Px = zeros(length(x), params.N);
Py = zeros(length(x), params.N);
for i = 1:length(x)
    if strcmp(basis, 'fourier')
        Px(i,:) = fourierLift( x(i, 1:n)', x(i, n+1:n+p)' )';
        Py(i,:) = fourierLift( y(i, 1:n)', y(i, n+1:n+p)' )';
    elseif strcmp(basis, 'poly')
        Px(i,:) = polyLift( x(i, 1:n)', x(i, n+1:n+p)' )';
        Py(i,:) = polyLift( y(i, 1:n)', y(i, n+1:n+p)' )';
    end
end

K = size(Px,1);
N = size(Px,2);

%% Solve for inital Koopman Operator wish subset of data points

% Build Apx sparsely with 10% of snapshotPairs
Ktithe = floor(K/10);   % roughly 10% of total snapshotPairs 
row = zeros(Ktithe*N^2,1);
col = zeros(Ktithe*N^2,1);
val = zeros(Ktithe*N^2,1);
for i = 1 : Ktithe*N
    row(1+(i-1)*N : i*N) = i;
    col(1+(i-1)*N : i*N) = mod(i-1,N) * N +1 : mod(i-1,N) * N + N;
    val((1+(i-1)*N : i*N)) = Px( ceil( i/N ) , : );
end
Apx = sparse(row, col, val, K*N, N*N);
bpy = reshape(Py', [K*N,1]);

% Call function that solves LP problem (EDIT HERE!!!)
% (Call Ram's function, something like this)
Uvec = robustKoopmanLP(sparse( Apx ), ( bpy ) );

%% Check which points the solution holds for, and repeat process as necessary

optimal = false;
Kcurrent = Ktithe;
while ~optimal
    for i = 1:K
        if ~any(find(satisfied == i))   % only check points that weren't part of the optimization problem
            
            % check if the constraints are satisfied for this point
            
        end
    end
end

%% Reshape the Koopman operator as a vector
U = reshape(Uvec, [N,N]);



end

