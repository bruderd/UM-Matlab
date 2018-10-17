function U = get_KoopmanConstGen( x,y, params )
%get_KoopmanConstGen: Find the best possible koopman operator given
%snapshot pairs using constraint generation to deal with large data sets.
%   Detailed explanation goes here

[n, p] = deal(params.n, params.p);

Px = zeros(length(x), params.N);
Py = zeros(length(x), params.N);
for i = 1:length(x)
    Px(i,:) = stateLift( x(i, 1:n)', x(i, n+1:n+p)' )';
    Py(i,:) = stateLift( y(i, 1:n)', y(i, n+1:n+p)' )';
end

K = size(Px,1);
N = size(Px,2);

%% Solve for inital Koopman Operator wish subset of data points

% Build Apx sparsely with 10% of snapshotPairs
Ktithe = min( floor(K/2) , 1000 );   % roughly 50% of total snapshotPairs, at most 1000 
row = zeros(Ktithe*N^2,1);
col = zeros(Ktithe*N^2,1);
val = zeros(Ktithe*N^2,1);
for i = 1 : Ktithe*N
    row(1+(i-1)*N : i*N) = i;
    col(1+(i-1)*N : i*N) = mod(i-1,N) * N +1 : mod(i-1,N) * N + N;
    val((1+(i-1)*N : i*N)) = Px( ceil( i/N ) , : );
end
Apx = sparse(row, col, val, Ktithe*N, N*N);
bpy = reshape(Py(1:Ktithe,:)', [Ktithe*N,1]);

% Call function that solves QP problem
Uvec = solve_KoopmanQP(Apx, bpy, params);
U = reshape(Uvec, [N,N]);

%% Check which points the solution holds for, and repeat process as necessary

optimal = false;
satConst = zeros(K,1);      % logical array that stores which constraints are satisfied
while ~optimal
    unsatisfied = [];       % stores the indices of all the unsatisfied constraints
    count = 0;
    
    % check if constraints are satisfied
    for i = 1:K
        if satConst(i,1) == 0         %~any(find(satisfied == i))   % only check points that weren't part of the last optimization problem
            satConst(i,1) = all( Px(i,:) * U - Py(i,:) <= params.epsilon');        %check if the constraints are satisfied for this point
            if satConst(i,1) == 0
                unsatisfied = [unsatisfied, i];     % store index of the unsatisfactory point
                count = count + 1;
            end
        end
        % ensures we add at most 1000 new points to our optimization problem
        if count > 1000
            break;
        end
    end
    
    % print progress and check if all (most) of the constraints are satisfied
    progress = 100 * sum(satConst)/K;
    disp(['Epsilon = ', num2str(mean(params.epsilon))]);     % print the average value of epsilon
    disp([num2str(progress), '% of constraints satisfied']);    % print percentage of constraints satisfied
    if (sum(satConst) > 0.9*K)
        optimal = true;
        break;
    end
    
    % Construct new A and b matrices
    Kunsat = length(unsatisfied);   % number of points where constraints unsatisfied
    row = zeros(Kunsat*N^2,1);
    col = zeros(Kunsat*N^2,1);
    val = zeros(Kunsat*N^2,1);
    for i = 1 : Kunsat
        for j = 1 : N
            index = N*(i-1) + j;
            row(1+(index-1)*N : index*N) = index;
            col(1+(index-1)*N : index*N) = mod(index-1,N) * N + 1 : mod(index-1,N) * N + N;
            val((1+(index-1)*N : index*N)) = Px( unsatisfied(i) , : );
        end
    end
    Apx_add = sparse(row, col, val, Kunsat*N, N*N);
    bpy_add = reshape( Py(unsatisfied,:)' , [Kunsat*N,1]);
    Apx = [Apx; Apx_add];
    bpy = [bpy; bpy_add];
    
    % Call function to solve QP
    Uvec = solve_KoopmanQP(Apx, bpy, params);
    U = reshape(Uvec, [N,N]); 
    
end


end

