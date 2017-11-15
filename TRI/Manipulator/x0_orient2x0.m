function x0 = x0_orient2x0( x0_orient, params )
%x0_orient2x: Converts the orientation in global coordinates to position
%and orientation in local coordinates
%   Detailed explanation goes here

p = params.p;       % total number of modules

%% calculate x in terms of x0_orient
x = zeros(6*p, 1);

x(4:6, 1) = x0_orient(1:3,1);
x(1:3, 1) = euler2cart(x0_orient(1:3,1), params.L(1));
for i = 2:p
   x0i_orient = x0_orient(1+3*(i-1) : 3*i, 1);     % state of ith module in global coordinates
   x0im1_orient = x0_orient(1+3*(i-2) : 3*(i-1), 1);     % state of (i-1)th module in global coordinates
   
   R0_i = Reuler(x0i_orient);       % rotation matrix from ith to 0 frame
   R0_im1 = Reuler(x0im1_orient);   % rotation matrix from (i-1)th to 0 frame
%    Rim1_0 = pinv(R0_im1);           % rotation matrix from to (i-1)th frame
   Rim1_0 = R0_im1';           % rotation matrix from to (i-1)th frame

   
   xi_orient = rot2euler_sym(Rim1_0 * R0_i);
   xi(1:3, 1) = euler2cart(xi_orient, params.L(i));      % set position component of xi
   xi(4:6, 1) = xi_orient;   % set orientation component of xi
   
   x(1+6*(i-1) : 6*i, 1) = xi;
end

x0 = x2x0(x, params);

end

