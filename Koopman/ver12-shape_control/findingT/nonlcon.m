function [ c , ceq ] = nonlcon( arg , params)
%nonlcon: Nonlinear constraints
%   arg - optimization decision variable
%   params - system parameters

N = params.N;

% convert decision variable to T matrix
T = reshape( arg( 1 : N^2 ) , [N,N] );

% Inequality constraints ( c <= 0 )
c = [ 1e-3 - abs( det( T ) );...    % det(T) ~= 0
     abs( det(T) ) - 2 ];          % |det(T)| <= 2

% Equality constraints ( ceq = 0 )

% Set the first row equal to known Taylor series coefficients
ceq = [];

end

