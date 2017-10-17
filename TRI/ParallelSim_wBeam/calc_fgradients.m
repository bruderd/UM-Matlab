function [ output_args ] = calc_fgradients(  )
%calc_fgradients
%   WORK IN PROGRESS, WILL NEED LATER FOR OPTIMIZATION

syms psi theta phi dpsi dtheta dphi ddpsi ddtheta ddphi L n a1 a2 a3 a4 b1 b2 b3 b4 N1 N2 N3 N4 

x = [psi, theta, phi, dpsi, dtheta, dphi]';
xdot = [dpsi, dtheta, dphi, ddpsi, ddtheta, ddphi]';

f = [xdot(1) - x(4);...
     xdot(2) - x(5);...
     xdot(3) - x(6);...
     tau(1) - 
     tau(2) - 
     tau(3) - 
     ];
end

