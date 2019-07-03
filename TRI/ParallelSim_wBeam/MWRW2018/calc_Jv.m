function Jv = calc_Jv( x, params )
%calc_Jv
%   Defines the Jacobian, Jv, such that fq = Jv * u

n = params.numFREEs;   % the number of FREE actuators in system
B = params.B;   % FREE fiber length (array)
N = params.Nf;   % FREE total fiber windings in revolutions (when relaxed) (array)
L = params.Lspine;   % lenght of central spine

% Define the extension and twist of each FREE with respect to x
q = euler2free(x,params);
s = q(1:n);
w = q(n+1:2*n);

Jv = zeros(2*n,n);      % initialize Jv
for i = 1:n
    Jv(i,i) = (pi*(B(i)^2 - 3*(L+s(i))^2) / (2*pi*N(i)+w(i))^2);   
    Jv(n+i,i) = 2*pi*(L+s(i))*((L+s(i))^2 - B(i)^2) / (2*pi*N(i)+w(i))^3;
end

end

