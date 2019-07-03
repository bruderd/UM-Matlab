function xcart = euler2cart( xeul, params )
%Gives the cartesian coordinates of the end effector induced by the orientation under
%assumption 2.1
%   Detailed explanation goes here

[psi,theta,phi] = deal(xeul(1), xeul(2), xeul(3));
L = params.Lspine;

xcart(1:3,1) = [(L/2)*(cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi));...
               (L/2)*(sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi));...
               (L/2)*(cos(theta)*cos(psi)-1)];

% NEED TO WRITE CALC_JX FUNCTION FOR THIS PART TO WORK
% Jx = calc_Jx(xeul, params);
% 
% xcart(4:6,1) = Jx * xeul(4:6);

end

