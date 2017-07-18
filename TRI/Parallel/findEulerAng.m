function eulerAng = findEulerAng( P, params )
%Finds the Euler angles describing the orientation of end effector given
%pressure in each FREE. This assumes there is no external load applied.
%   Detailed explanation goes here
          
eulerAng = lsqnonlin(@(eulerAng) Meff(eulerAng, P, params), [0;0;0]);

% calculate the remainder of the optimization (really just for debugging)
remainder = Meff(eulerAng, P, params);

end

