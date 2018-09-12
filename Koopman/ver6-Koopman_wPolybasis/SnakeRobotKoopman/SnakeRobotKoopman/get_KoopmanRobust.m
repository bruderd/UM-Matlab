function U = get_KoopmanRobust( x,y, params )
%get_KoopmanRobust: Uses robust linear regression to solve for Koopman operator
%   Detailed explanation goes here

[n, p] = deal(params.n, params.p);

Px = zeros(length(x), params.N);
Py = zeros(length(x), params.N);
for i = 1:length(x)
    Px(i,:) = polyLift( x(i, 1:n)', x(i, n+1:n+p)' )';
    Py(i,:) = polyLift( y(i, 1:n)', y(i, n+1:n+p)' )';
%     Px(i,:) = sinLift( x(i, 1:n)', x(i, n+1:n+p)' )';
%     Py(i,:) = sinLift( y(i, 1:n)', y(i, n+1:n+p)' )';
end

% Vectorize/Matricizw Px and Py to pass to Ram's function
K = size(Px,1);
N = size(Px,2);

% Apx = zeros(K*N, N*N);
% for i = 1 : size(Px,1)
%     A = Px(i,:);
%     ACell = repmat({A}, 1, N);
%     BigA = blkdiag(ACell{:});
%     
%     Apx(1 + (i-1)*N : i*N, :) = BigA; 
% end

% Build Apx sparsely
row = zeros(K*N^2,1);
col = zeros(K*N^2,1);
val = zeros(K*N^2,1);
for i = 1 : K*N
    row(1+(i-1)*N : i*N) = i;
    col(1+(i-1)*N : i*N) = mod(i-1,N) * N +1 : mod(i-1,N) * N + N;
    val((1+(i-1)*N : i*N)) = Px( ceil( i/N ) , : );
end
Apx = sparse(row, col, val, K*N, N*N);

bpy = reshape(Py', [K*N,1]);

%% Call Ram's function
Uvec = robustKoopmanLP(sparse( Apx ), ( bpy ) );

%% Reshape the Koopman operator as a vector
U = reshape(Uvec, [N,N]);


end

