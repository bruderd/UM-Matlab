function Meff = Meff_v2( eulerAng, P, params )
% Meff_v2.m
%   Determine the net torques acting on the end effector. This differs from
%   Meff in that it uses the full manipulator Jacobian to calculate the
%   torques.

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
    s(i) = sqrt((L - sin(theta)*a(i) + cos(theta)*sin(psi)*b(i))^2 + (a(i)^2 + b(i)^2)*phi^2) - L;    % technically this only holds if FREE initial lengths are the same as spine.
end
  
% Calculate the jacobian that depends on the values of s and w
Jdagger(1,:) = [(a(1)+b(1))^2*phi/(s(1)+L1),...
                (L + b(1)*cos(theta)*sin(psi) - a(1)*sin(theta))*(-a(1)*cos(theta) - b(1)*sin(psi)*sin(theta))/(s(1) + L1),...
                b(1)*cos(psi)*cos(theta)*(L + b(1)*cos(theta)*sin(psi) - a(1)*sin(theta))/(s(1) + L1)];
Jdagger(2,:) = [(a(2)+b(2))^2*phi/(s(2)+L2),...
                (L + b(2)*cos(theta)*sin(psi) - a(2)*sin(theta))*(-a(2)*cos(theta) - b(2)*sin(psi)*sin(theta))/(s(2) + L2),...
                b(2)*cos(psi)*cos(theta)*(L + b(2)*cos(theta)*sin(psi) - a(2)*sin(theta))/(s(2) + L2)];
Jdagger(3,:) = [(a(3)+b(3))^2*phi/(s(3)+L3),...
                (L + b(3)*cos(theta)*sin(psi) - a(3)*sin(theta))*(-a(3)*cos(theta) - b(3)*sin(psi)*sin(theta))/(s(3) + L3),...
                b(3)*cos(psi)*cos(theta)*(L + b(3)*cos(theta)*sin(psi) - a(3)*sin(theta))/(s(3) + L3)];            
Jdagger(4,:) = [(a(4)+b(4))^2*phi/(s(4)+L4),...
                (L + b(4)*cos(theta)*sin(psi) - a(4)*sin(theta))*(-a(4)*cos(theta) - b(4)*sin(psi)*sin(theta))/(s(4) + L4),...
                b(4)*cos(psi)*cos(theta)*(L + b(4)*cos(theta)*sin(psi) - a(4)*sin(theta))/(s(4) + L4)];
Jdagger(5:8,:) = [1 0 0; 1 0 0; 1 0 0; 1 0 0];   

partials = [-pi*(B1^2 - 3*(L1+s(1))^2) / (2*pi*N1+w(1))^2, 0, 0, 0;...
            0, -pi*(B2^2 - 3*(L2+s(2))^2) / (2*pi*N2+w(2))^2, 0, 0;...
            0, 0, -pi*(B3^2 - 3*(L3+s(3))^2) / (2*pi*N3+w(3))^2, 0;...
            0, 0, 0, -pi*(B4^2 - 3*(L4+s(4))^2) / (2*pi*N4+w(4))^2;...
            2*pi*(L1+s(1))*((L1+s(1))^2 - B1^2) / (2*pi*N1+w(1))^3, 0, 0, 0;...
            0, 2*pi*(L2+s(2))*((L2+s(2))^2 - B2^2) / (2*pi*N2+w(2))^3, 0, 0;...
            0, 0, 2*pi*(L3+s(3))*((L3+s(3))^2 - B3^2) / (2*pi*N3+w(3))^3, 0;...
            0, 0, 0, 2*pi*(L4+s(4))*((L4+s(4))^2 - B4^2) / (2*pi*N4+w(4))^3];
            
            
Meff = Jdagger'*partials*P + Ks*(-s) + Kw*w;


end