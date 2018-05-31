function Jq = calcJq( q, params )
%calcJqi: Calculates the fluid Jacobian for a single FREE
%   Detailed explanation goes here

num = params.num;
[L, B, N] = deal(params.L, params.B, params.N);

dl = q(1:2:end);
dphi = q(2:2:end);

% Define volume Jacobain
Jq = zeros(sum(num), 2*sum(num));      % initialize Jv
for i = 1:num
    Jq_i = [(pi*(B(i)^2 - 3*(L(i)+dl(i))^2) / (2*pi*N(i)+dphi(i))^2),...
             2*pi*(L(i)+dl(i))*((L(i)+dl(i))^2 - B(i)^2) / (2*pi*N(i)+dphi(i))^3];
         
    Jq( i, 2*(i-1)+1 : 2*i ) = Jq_i;     % stack volume jacobian for each actuator diagonally to form Jv.   
end

end

