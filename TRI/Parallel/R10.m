function R = R10( phi, theta, psi )
%Calculates the rotation matrix to transform form coordinates in frame-1 to
%frame-0, given the Euler angles describing their relative orientation.
%   Detailed explanation goes here

R = [cos(phi)*cos(theta), cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi), cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi);...
     sin(phi)*cos(theta), sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi), sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi);...
     -sin(theta), cos(theta)*sin(psi), cos(theta)*cos(psi)];

end

