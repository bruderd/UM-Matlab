function [ JxT ] = getJx_symbolic( num, a, d, Gama, R, L )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% 
% Example Actuator parameters:
% num = 3;    % number of FREEs in combination
% Gama = deg2rad([15, -15, 52.77]); % relaxed fiber angle of each FREE
% R = (12.319e-3)/2 * ones(1,num);  % relaxed radius of each FREE [m]
% L = 0.16 * ones(1,num);   %  relaxed length of each FREE [m] 
% d = [0, -0.02286, 0 ; 0.01979676, 0.01143, 0 ; -0.01979676, 0.01143, 0]'; % location of attachment points to the end effector [m]
% a = [0,0,1 ; 0,0,1 ; 0,0,1]';    % direction of FREE axis at attachment point [unit vector]

B = abs(L ./ cos(Gama));   % fiber length (must be positive))
N = -L ./ (2*pi*R) .* tan(Gama); % total fiber windings in revolutions (when relaxed)

q = sym('q', [2*num, 1], 'real');    % symbolic FREE variable
x = sym('x', [6,1], 'real');    % symbolic state variable

% 2 DOF rotation and translation rig (axes not aligned)
for i = 1 : num
    q(2*(i-1) + 1) = -L(i) + sqrt( (L(i)+x(3))^2 + ( x(4)*(d(1,i) + d(2,i)) )^2);
    q(2*(i-1) + 2) = x(4);
end

% matrix representation for cross product with ith column of d 
dcross = zeros(3,3,num);
for i = 1 : num
    dcross(:,:,i) = [0, -d(3,i), d(2,i);...
                     d(3,i), 0, -d(1,i);...
                     -d(2,i), d(1,i), 0];
end

% transformation matrix for each FREE
D = zeros(6,2,num);
for i = 1 : num
   D(:,:,i) = [...
                [a(:,i) , zeros(3,1)];...
                [dcross(:,:,i) * a(:,i) , zeros(3,1)] + [zeros(3,1) , a(:,i)]
              ];
end

% fluid jacobian for each FREE
Jq = sym('Jq', [1,2,num], 'real');
for i = 1 : num
    q1 = q(2*(i-1) + 1);
    q2 = q(2*(i-1) + 2);
    Jq(1,1,i) = pi * ( B(i)^2 - 3*(L(i)+q1)^2 ) / (2 * pi * N(i) + q2)^2;
    Jq(1,2,i) = 2 * pi * (L(i)+q1) * ((L(i)+q1)^2 - B(i)^2) / (2 * pi * N(i) + q2)^3;
end

% Jx for each FREE
Jxi = sym('Jxi', [1,6,num], 'real');
for i = 1 : num
    Jxi(:,:,i) = Jq(:,:,i) * D(:,:,i)' ;
end

% Combine all Jxi's together into a single jacobian Jx
JxT = Jxi(:,:,1)';
for i = 2:num
    JxT = [JxT, Jxi(:,:,i)'];
end

end

