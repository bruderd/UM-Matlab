function [ c , ceq ] = nonlcon( arg , params)
%nonlcon: Nonlinear constraints
%   arg - optimization decision variable
%   params - system parameters

N = params.N;

% convert decision variable to T matrix
T = reshape( arg( 1 : N^2 ) , [N,N] );

% Inequality constraints ( c <= 0 )
c = [ 1e-3 - abs( det( T ) );...    % det(T) ~= 0
     abs( det(T) ) - 10 ];          % |det(T)| <= 10

% Equality constraints ( ceq = 0 )
ceq = [];

% % Set the first row equal to known Taylor series coefficients (specific to
% % poly2 basis functions)
% ceq = [ arg(1) + 3 ;...
%         arg(7);...
%         arg(13);...
%         arg(19) - 2;...
%         arg(25);...
%         arg(31)];

% % set determinant equal to 1
% ceq = det(T) - 1;

end

