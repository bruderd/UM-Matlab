function [ c , ceq ] = nonlcon_transmodel( arg , params)
%nonlcon: Nonlinear constraints
%   arg - optimization decision variable
%   params - system parameters

N = params.N;
n = params.n;

% convert decision variable to T matrix
T = reshape( arg , [n,N] );

% Inequality constraints ( c <= 0 )
% c = [ 1e-3 - abs( det( T * T' ) );...    % det(T) ~= 0
%       abs( det( T * T' ) ) - 10 ;...          % |det(T)| <= 10
%       abs( mean(mean( abs(T(:,1:end-1)) ) ) - mean( abs(T(:,end)) ) ) - 1e-1];    % last column should be same order of magnitude as other columns 
% c = [ 1e-3 - abs( det( T * T' ) );...    % det(T) ~= 0
%       1e-3 - norm(T);...
%       norm(T) - 1e0;...
%       abs( mean(mean( abs(T(:,1:end-1)) ) ) - mean( abs(T(:,end)) ) ) - 1e-1];    % last column should be same order of magnitude as other columns 
c = [ 1e-3 - abs( det( T * T' ) );...    % det(T) ~= 0
      abs( det( T * T' ) ) - 10 ];          % |det(T)| <= 10
c = [];
    
cor_mtx = T*T';
cor_vec = cor_mtx(:);
eye_mtx = eye(n);
eye_vec = eye_mtx(:);
  
% Equality constraints ( ceq = 0 )
ceq = [ T(:,end) ;... % last column should be zeros
        cor_vec - eye_vec];

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