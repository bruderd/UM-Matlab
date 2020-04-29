function [ yplus , xplus ] = vf_bilinear(k,x,u,A,N,B,C)
%vf_bilinear: Dynamics of sample discrete dynamical system (SISO)
%   Detailed explanation goes here

xplus = A * x + N * x * u(k) + B * u(k); 
yplus = C * xplus;   % nonlinear output map

end