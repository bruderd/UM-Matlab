function euler = rot2euler( R )
%rot2euler: Exracts intrinsic ZYX euler angles from a 3x3 rotation matrix
%   Detailed explanation goes here

euler = zeros(3,1);

if R(1,1) == 0 && R(2,1) == 0
    euler(1) = atan2(R(1,2), R(2,2));
    euler(2) = pi/2;
    euler(3) = 0;
else
    euler(1) = atan2(R(3,2), R(3,3));
    euler(2) = atan2(-R(3,1), sqrt(R(1,1)^2 + R(2,1)^2));
    euler(3) = atan2(R(2,1), R(1,1));
end

end

