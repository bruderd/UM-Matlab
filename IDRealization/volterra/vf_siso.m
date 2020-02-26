function [ yplus , xplus ] = vf_siso(k,x,u)
%vf_siso: Dynamics of sample discrete dynamical system (SISO)
%   Detailed explanation goes here


A = [ 1 , 0.01 ; -0.01 , 1]; B = [0 ; 0.01];

xplus = A * x + B * u(k);  % linear dynamics

yplus = sin( xplus(1) ) + xplus(1)*xplus(2);   % nonlinear output map

end

