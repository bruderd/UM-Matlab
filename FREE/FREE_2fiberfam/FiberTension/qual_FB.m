% qual_FB.m
%   Uses the force balance equations to qualitatively determine the initial
%   behavior of a FREE when pressurized
%
%   Calls solveFB.m to determine values of parameters at pressure P_test
%   Or you can change it to call solveFB_v2.m to solve newer version of it.

function [dL, dphi] = qual_FB(P_test, x_rest)

% Give names to x_rest parameters
[P0, gama0, betta0, r0, L0, phi0] = deal(x_rest(1), x_rest(2), x_rest(3), x_rest(4), x_rest(5), x_rest(6));

% Solve for Tgama, Tbetta, P, gama, betta, r, L, phi
[T_gama, T_betta, P, gama, betta, r, L, phi] = solveFB_v2(P_test, x_rest);

% Determine the direction of deformation:

% % Quick test: What if I fix the tensions for all anlges?
% P = P_test;
% T_gama = 0;
% T_betta = 2*P*pi*r0^2/(4*cos(deg2rad(54.7)));

% % via forces:
% dL = sign(P*pi*r0^2 - 2*(T_gama*cos(gama0) + T_betta*cos(betta0)));
% 
% if abs(phi) < 1e-8
%     dphi = 0;
% else
%     dphi = sign(2*r0*(T_gama*sin(gama0) + T_betta*sin(betta0)));
% end

% via displacements:
dL = sign(L - x_rest(5));
if abs(phi) < 1e-8
    dphi = 0;
else
    dphi = sign(phi);
end

end

