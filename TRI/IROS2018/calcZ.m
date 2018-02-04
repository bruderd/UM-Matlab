function Z = calcZ( q, P, params )
%maxZ: Calculates the maximum Z = [F,M]' for each free
%   Zmax is a vertical concatenation or Zmax_i for each FREE.
%   Pmax should be a vector of dim n.

num = params.num;   % number of FREEs in parallel
B = params.B;   % FREE fiber length (array)
N = params.N;   % FREE total fiber windings in revolutions (when relaxed) (array)
L = params.L;   % length of each FREE

s = q(1:2:end);
w = q(2:2:end);

% Define volume Jacobain
Jv = zeros(sum(num), 2*sum(num));      % initialize Jv
for i = 1:num
    Jv_ki = [(pi*(B(i)^2 - 3*(L(i)+s(i))^2) / (2*pi*N(i)+w(i))^2),...
             2*pi*(L(i)+s(i))*((L(i)+s(i))^2 - B(i)^2) / (2*pi*N(i)+w(i))^3];
         
    Jv( i, 2*(i-1)+1 : 2*i ) = Jv_ki;     % stack volume jacobian for each actuator diagonally to form Jv.   
end

Zmax = Jv' * P;

end