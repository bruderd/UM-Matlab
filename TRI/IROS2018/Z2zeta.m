function zeta = Z2zeta( Z, params )
%Z2zeta: Converts FREE forces in FREE coordinates to end effector
%coordinates
%   Detailed explanation goes here

num = params.num;   % number of FREEs in parallel
d = params.d;
p = params.p;

for i = 1:num
    diX = [0 , -d(3,i), d(2,i); d(3,i), 0, -d(1,i); -d(2,i), d(1,i), 0];
    Di = [p(:,i), zeros(3,1); [zeros(3,1), p(:,i)] + diX*[p(:,i), zeros(3,1)]];
    zeta(6*(i-1)+1 : 6*i,1) = Di * Z(2*(i-1)+1 : 2*i, 1);
end

end

