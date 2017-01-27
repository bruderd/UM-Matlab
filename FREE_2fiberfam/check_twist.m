%% Checks that the twist at end due to one fiber is the same as twist due to another fiber
%       This can only be run after ode15i of a 2fiber system

gama0 = x0(2);
betta0 = x0(3);
r0 = x0(4);
L0 = x0(5);

gama = y(end,2);
betta = y(end,3);
r = y(end,4);
L = y(end,5);

phi_gama = (tan(gama)*L/r - tan(gama0)*L0/r0);
phi_betta = (tan(betta)*L/r - tan(betta0)*L0/r0);