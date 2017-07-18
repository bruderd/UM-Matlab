function Meff = Meff( eulerAng, P, params )
% netF.m
%   Determine the net force and torque on the central element 

% rename relevant parameters

[phi, theta, psi] = deal(eulerAng(1), eulerAng(2), eulerAng(3));

R = R10(phi, theta, psi);

L = params.Lspine;
[L1, B1, N1] = deal(params.L1, params.B1, params.N1); 
[L2, B2, N2] = deal(params.L2, params.B2, params.N2);
[L3, B3, N3] = deal(params.L3, params.B3, params.N3); 
[L4, B4, N4] = deal(params.L4, params.B4, params.N4);

a = [params.attach1(1), params.attach2(1), params.attach3(1), params.attach4(1)];
b = [params.attach1(2), params.attach2(2), params.attach3(2), params.attach4(2)];

Ks = [params.kelast1(1) * b(1), params.kelast2(1) * b(2), params.kelast3(1) * b(3), params.kelast4(1) * b(4);...
      params.kelast1(1) * a(1), params.kelast2(1) * a(2), params.kelast3(1) * a(3), params.kelast4(1) * a(4);...
      0, 0, 0, 0];

Kw = [0, 0, 0, 0;...
      0, 0, 0, 0;...
      params.kelast1(2), params.kelast2(2), params.kelast3(2), params.kelast4(2)];

% Calculate s and w from Euler angles
w = phi * ones(params.numFREEs, 1);
s = zeros(params.numFREEs, 1);
for i = 1:params.numFREEs
    s(i) = norm( R*[a(i); b(i); 0] - [a(i); b(i); - L] ) - L;
%     s(i) = sqrt((L - sin(theta)*a(i) + cos(theta)*sin(psi)*b(i))^2 + (a(i)^2 + b(i)^2)*phi^2) - L;    % technically this only holds if FREE initial lengths are the same as spine.
end
  
% Calculate the jacobian that depends on the values of s and w
J(1,:) = [-pi*(B1^2 - 3*(L1+s(1))^2) / (2*pi*N1+w(1))^2 * b(1),...
          -pi*(B2^2 - 3*(L2+s(2))^2) / (2*pi*N2+w(2))^2 * b(2),...
          -pi*(B3^2 - 3*(L3+s(3))^2) / (2*pi*N3+w(3))^2 * b(3),...
          -pi*(B4^2 - 3*(L4+s(4))^2) / (2*pi*N4+w(4))^2 * b(4)];

J(2,:) = [-pi*(B1^2 - 3*(L1+s(1))^2) / (2*pi*N1+w(1))^2 * a(1),...
          -pi*(B2^2 - 3*(L2+s(2))^2) / (2*pi*N2+w(2))^2 * a(2),...
          -pi*(B3^2 - 3*(L3+s(3))^2) / (2*pi*N3+w(3))^2 * a(3),...
          -pi*(B4^2 - 3*(L4+s(4))^2) / (2*pi*N4+w(4))^2 * a(4)];
      
J(3,:) = [2*pi*(L1+s(1))*((L1+s(1))^2 - B1^2) / (2*pi*N1+w(1))^3,...
          2*pi*(L2+s(2))*((L2+s(2))^2 - B2^2) / (2*pi*N2+w(2))^3,...
          2*pi*(L3+s(3))*((L3+s(3))^2 - B3^2) / (2*pi*N3+w(3))^3,...
          2*pi*(L4+s(4))*((L4+s(4))^2 - B4^2) / (2*pi*N4+w(4))^3];

      

Meff = J*P + Ks*(-s) + Kw*w;


end