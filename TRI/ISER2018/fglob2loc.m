function floc = fglob2loc( fglob, x, params )
%fglob2loc: Expresses forces in the local end effector body fixed
%coordinates that were originally expressed in global space-fixed world
%coordinates.
%   This function is especially useful for gravitational forces

% convert euler angles into a rotation matrix
R = eul2rotm(x(4:6));
Rtot = blkdiag(R,R);

floc = Rtot * fglob;

end

