function [ Uvec , epsilon ] = robustKoopmanLP_withCG( A , b , params )
%robustKoopmanLP_withCG: Summary of this function goes here
%   Detailed explanation goes here

na = size( A, 1 );
nx = size( A, 2 );
nb = size( b, 1 );

Istacked = repmat( speye(params.N) , (size( A,1 ) / params.N) , 1);

tildeA = [ A , -Istacked, zeros(na,1) ;...
          -A , -Istacked, zeros(na,1) ;...
          ones(1,params.N^2) , zeros(1,params.N), -1 ;...
          -ones(1,params.N^2) , zeros(1,params.N), -1];
tildeb = [b ; -b; 0 ; 0];

f = zeros( nx + params.N + 1, 1 );
f( (nx + 1):(nx + params.N) ) = 1;
f( end ) = 100;     % penalty for size of elements in U

[ xout, fval, exitflag ] = linprog(f, tildeA, tildeb);

Uvec = xout(1:nx);
epsilon = xout( (nx + 1):(nx + params.N) );

end

