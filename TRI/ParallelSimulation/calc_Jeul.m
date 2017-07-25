function Jeul = calc_Jeul( x, params )
%cacl_Jeul
%   Calulate the Jacobian, Jeul, such that: qdot = Jeul * [dpsi,dtheta,dphi]'

n = params.numFREEs;   % the number of FREE actuators in system
a = params.xattach;   % x-cordinate of FREE attachment points (array)
b = params.yattach;   % y-coordinate of FREE attachment points (array)
L = params.Lspine;   % lenght of central spine

[psi, theta, phi] = deal(x(1), x(2), x(3));

Jeul = zeros(2*n,3);    % initialize matrix
for i = 1:n
   denom = sqrt((L - sin(theta)*a(i) + cos(theta)*sin(psi)*b(i))^2 + (a(i)^2 + b(i)^2)*phi^2);
   Jeul(i,1) = b(i)*cos(psi)*cos(theta)*(L + b(i)*cos(theta)*sin(phi) - a(i)*sin(theta)) / denom;
   Jeul(i,2) = (L + b(i)*cos(theta)*sin(psi) - a(i)*sin(theta))*(-a(i)*cos(theta) - b(i)*sin(psi)*sin(theta)) / denom;
   Jeul(i,3) = (a(i)^2 + b(i)^2)*phi / denom;
   Jeul(n+i, 3) = 1;
end

end

