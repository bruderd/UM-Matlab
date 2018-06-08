function q = euler2free_MWRW( x, params )
%euler2free
%   Converts the euler angles to a vector of extensions and twists of
%   FREEs. This works for an arbitrary number of FREE actuators.
%   NOTE: This was changed from the original version to be consistent with
%   the QP optimization code.

n = params.numFREEs;   % the number of FREE actuators in system
a = params.xattach;   % x-cordinate of FREE attachment points (array)
b = params.yattach;   % y-coordinate of FREE attachment points (array)
L = params.Lspine;   % lenght of central spine

[psi, theta, phi] = deal(x(1), x(2), x(3));

q = zeros(2*n,1);   % initialize q
for i = 1:n
   q(2*i-1) = -L + sqrt((L - sin(theta)*a(i) + cos(theta)*sin(psi)*b(i))^2 + (a(i)^2 + b(i)^2)*phi^2);
   q(2*i) = phi;
end

end

