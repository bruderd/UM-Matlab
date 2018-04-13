function y=rms(x,dim)
%Root Mean Square.
% For vectors, RMS(X) is the root mean square value of the elements in X.
%
% For matrices, RMS(X) is a row vector containing the root mean square
% value of each column.  For N-D arrays, RMS(X) is the root mean square
% value of the elements along the first non-singleton dimension of X.
%
% RMS(X,DIM) takes the root mean square along the dimension DIM of X. 
%
% See also MEAN, MEDIAN, STD, MIN, MAX, VAR, COV, MODE.

% D.C. Hanselman, University of Maine, Orono, ME 04469
% MasteringMatlab@yahoo.com
% Mastering MATLAB 7
% 2006-01-11

if nargin==1 && isvector(x)             % rms(x), speed for vector case
   y=norm(x)/sqrt(length(x));
elseif nargin==1                        % rms(x), general case
   dim=min(find(size(x)~=1));
   if isempty(dim)
      dim = 1;
   end
   y=sqrt(sum(x.^2)/size(x,dim));
else                                    % rms(x,dim), dimension specified
   y=sqrt(sum(x.^2,dim)/size(x,dim));
end