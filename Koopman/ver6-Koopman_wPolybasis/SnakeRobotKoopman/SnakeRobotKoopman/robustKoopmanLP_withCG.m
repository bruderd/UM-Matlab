function [ Uvec , epsilon ] = robustKoopmanLP_withCG( A , b , params )
%robustKoopmanLP_withCG: Summary of this function goes here
%   Detailed explanation goes here

na = size( A, 1 );
nx = size( A, 2 );
nb = size( b, 1 );

%% Ridge Method

% Istacked = repmat( speye(params.N) , (size( A,1 ) / params.N) , 1);
% 
% tildeA = [ A , -Istacked, zeros(na,1) ;...
%           -A , -Istacked, zeros(na,1) ;...
%           ones(1,params.N^2) , zeros(1,params.N), -1 ;...
%           -ones(1,params.N^2) , zeros(1,params.N), -1];
% tildeb = [b ; -b; 0 ; 0];
% 
% f = zeros( nx + params.N + 1, 1 );
% f( (nx + 1):(nx + params.N) ) = 1;  % penalty error (was 1)
% f( end ) = 0;     % penalty for size of elements in U (was 100)
% 
% [ xout, fval, exitflag ] = linprog(f, tildeA, tildeb);


%% Lasso method
% x is vectorized Koopman operator, decomposed into positive and negative parts of each entry x = [u11+, ..., uNN+, u11-, ... , uNN-]';
% Uvec = M * x, where M subtracts the + and - parts of each entry: uij+ - uij-

M = [eye(params.N^2) , -eye(params.N^2)];

% L2 error as cost function
H = M' * A' * A * M;
f = -M' * A' * b;

% L1 regularization enforced as constraint
t = params.t;
Aq = [ -eye(2*params.N^2) ; ones(1 , 2*params.N^2) ];
bq = [ zeros(2*params.N^2 , 1) ; t ];

% Solve the quadratic program
[ x, fval, exitflag ] = quadprog(H, f, Aq, bq);

% Recover Uvec from the optimization variable
xout = M * x;


%% Set output

Uvec = xout(1:nx);

epsilon = 10;
% epsilon = 10.05e-2;

% epsilon = xout( (nx + 1):(nx + params.N) );

end

