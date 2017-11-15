function x0 = x2x0_sym( x, params )
%x02x: Converts symbolic state in global coordinates to symbolic local module coordinates
%   Detailed explanation goes here

p = params.p;

% selection matrix to isolate position component of module state
Spos = [eye(3), zeros(3,3)];
Seul = [zeros(3,3), eye(3)];

%% calculate x0 in terms of x
x0 = zeros(size(x));
x0 = sym(x0);

% define x0 in terms of x
% x0(1:6,1) = x(1:6,1);
x0(4:6, 1) = x(4:6,1);
x0(1:3, 1) = euler2cart(x(4:6,1), params.L(1));

R0 = zeros(3*p,3);
R0 = sym(R0);
R0(1:3, 1:3) = Reuler(x(4:6, 1));
for i = 2:p
    xi = x(1+6*(i-1) : 6*i, 1);     % local coordinates of ith module
    x0im1 = x0(1+6*(i-2) : 6*(i-1), 1);     % global coordinates of (i-1)th module
    
    Rim1_i = Reuler(xi(4:6));   % rotation matrix from ith to (i-1)th frame
    R0_im1 = R0(1+3*(i-2) : 3*(i-1), 1:3); % rotation matrix from (i-1)th frame to 0 (global) frame
    R0_i = R0_im1 * Rim1_i; % rotation matrix from ith frame to 0 (global) frame

    % just added, to get everything written in terms of euler angles...
    xi(1:3,1) = euler2cart(xi(4:6, 1), L(i));     % writing x position in terms of its angles.
    
    % perform coordinate transformation of xi to x0i
    x0i(1:3, 1) = x0im1(1:3,1) + R0_im1 * xi(1:3, 1);
    x0i(4:6, 1) = rot2euler_sym(R0_i);
    
    x0(1+6*(i-1) : 6*i, 1) = x0i;       % concatenate x0's of each module
    
    R0(1+3*(i-1) : 3*i, 1:3) = R0_i;    % storing the local to global rotations in the R0 matrix

end


end

