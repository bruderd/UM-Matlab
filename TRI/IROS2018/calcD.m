function D = calcD( d, a )
%Caculates coordinate transformation matrix D, given offset d and
%orientation a of the attachment points of FREE to the end effector

dx = [0, -d(3), d(2); d(3), 0, -d(1); -d(2), d(1), 0];

Dtop = [a, zeros(3,1)];
Dbottom = [dx*a, zeros(3,1)] + [zeros(3,1), a];
D = [Dtop; Dbottom];

end

