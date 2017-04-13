%% static_balance_ICRA
function [F] = static_balance_ICRA(P,x)
%Solve for all states given in input pressure
%   Detailed explanation goes here

% L0 = 3.152;   %short length FREE
% L0 = 7.25;    %long length FREE
L0 = 5.68;      %medium length FREE
% L0 = 5.808;     %70deg FREE
% L0 = 5.902;     %20deg FREE

gama0 = 0.7037; %40deg FREE
%gama0 = deg2rad(70);    %70deg FREE
%gama0 = deg2rad(20);    %20deg FREE

r0 = 3/16;

gama = x(1);
r = x(2);
L = x(3);
theta = L*tan(gama)/r;
phi = L*tan(gama)/r - L0*tan(gama0)/r0;

Fs = 1*(L);
Ms = 1*(theta);

F(1) = P*pi*r^2 - 2*P*pi*r^2*cot(gama)^2 + Fs;
F(2) = -2*P*pi*r^3*cot(gama) + Ms;        %no load
F(3) = cos(gama) - (L/L0)*cos(gama0);

