function coordinates = orient2pos( eulerAng, params )
%Gives the coordinates of the end effector induced by the orientation under
%assumption 2.1
%   Detailed explanation goes here

[phi,theta,psi] = deal(eulerAng(1), eulerAng(2), eulerAng(3));
L = params.Lspine;
R = R10(phi,theta,psi);

m = R*[0; 0; 1/L] - [0; 0; 1/L];    % slope of curve function
b = [0; 0; 1];                      % "initial value" of curve function

funx = @(l) m(1)*l + b(1);
funy = @(l) m(2)*l + b(2);
funz = @(l) m(3)*l + b(3);

x_eff = integral(funx,0,L);
y_eff = integral(funy,0,L);
z_eff = integral(funz,0,L);

coordinates = [x_eff; y_eff; z_eff] - [0; 0; L];

end

