function x = robustKoopmanLP( A, b )

nx = size( A, 2 );
nb = size( b, 1 );

% tildeA = sparse( 2 * nb + 2 * nx, nx + nb + 1 );
tildeb = zeros( 2 * nb + 2 * nx, 1 );

[ rowA, colA, sA ] = find( A );

rowtilA = rowA;
coltilA = colA;
stilA = sA;

rowtilA = [ rowtilA; nb + rowtilA ];
coltilA = [ coltilA; coltilA ];
stilA = [ stilA; -stilA ];

rowtilA = [ rowtilA; (1:nb)' ];
coltilA = [ coltilA; ( ( nx + 1 ):(nx + nb) )' ];
stilA = [ stilA; -1 * ones( nb, 1 ) ];

rowtilA = [ rowtilA; ( (nb + 1):2*nb )' ];
coltilA = [ coltilA;  ( ( nx + 1 ):(nx + nb) )' ];
stilA = [ stilA; -1 * ones( nb, 1 ) ];

rowtilA = [ rowtilA; ( 2*nb + ( ( 1:nx ) - 1 ) * 2 + 1 )' ];
coltilA = [ coltilA; ( ( ( 1:nx ) - 1 ) + 1 )' ];
stilA = [ stilA; 1 * ones( nx, 1 ) ];

rowtilA = [ rowtilA; ( 2*nb + ( ( 1:nx ) - 1 ) * 2 + 2 )' ];
coltilA = [ coltilA; ( ( ( 1:nx ) - 1 ) + 1 )' ];
stilA = [ stilA; -1 * ones( nx, 1 ) ];

rowtilA = [ rowtilA; ( 2*nb + ( ( 1:nx )- 1 ) * 2 + 1 )' ];
coltilA = [ coltilA; (nx + nb + 1) * ones( nx, 1 ) ];
stilA = [ stilA; -1 * ones( nx, 1 ) ];

rowtilA = [ rowtilA; ( 2*nb + ( ( 1:nx ) - 1 ) * 2 + 2 )' ];
coltilA = [ coltilA; (nx + nb + 1)  * ones( nx, 1 ) ];
stilA = [ stilA; -1 * ones( nx, 1 ) ];
%
% tildeA( 1:nb, 1:nx ) = A;
% tildeA( (nb + 1):2*nb, 1:nx ) = -A;
%
% for i = 1:nb
%     tildeA( i, nx + i ) = -1;
%     tildeA( nb + i, nx + i ) = -1;
% end

tildeb( 1:nb ) = b;
tildeb( (nb + 1):2*nb ) = -b;

% for i = 1:nx
%     tildeA( 2*nb + ( i - 1 ) * 2 + 1, ( i - 1 ) + 1 ) = 1;
%     tildeA( 2*nb + ( i - 1 ) * 2 + 2, ( i - 1 ) + 1 ) = -1;
%     tildeA( 2*nb + ( i - 1 ) * 2 + 1, end ) = -1;
%     tildeA( 2*nb + ( i - 1 ) * 2 + 2, end ) = -1;
% end

f = zeros( nx + nb + 1, 1 );
f( (nx + 1):(nx + nb) ) = 1;
f( end ) = 1e-6;

[ xout, fval, exitflag ] = linprog(f, sparse( rowtilA, coltilA, stilA, 2 * nb + 2 * nx, nx + nb + 1 ), tildeb);
% [ xout, fval, exitflag ] = linprog(f, tildeA, tildeb);


x = xout( 1:nx );
% sum( xout( nx+1:(nx + nb ) ) )
% xout( end )
% fval