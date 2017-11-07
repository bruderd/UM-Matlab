function R = Reuler( euler )
%Reuler: Calculates intrinsic ZYX rotation given Euler angles euler = [psi,theta, phi]'   
%   phi is angle about z-axis, theta is angle about y-axis, and psi is angle about x-axis.

if length(euler) ~= 3
    error('Input muse be a 3x1 vector of Euler angles');
end

[psi, theta, phi] = deal(euler(1), euler(2), euler(3));

R = [cos(phi)*cos(theta), cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi), cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi);...
     sin(phi)*cos(theta), sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi), sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi);...
     -sin(theta), cos(theta)*sin(psi), cos(theta)*cos(psi)];

end

