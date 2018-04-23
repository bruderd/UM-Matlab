function f = calcf( x, p, params )
%calcf: Calculates the force felt by the end effector with respect to its
%local body-fixed coordinate frame due to fiber forces only (no elastomer)
%   Detailed explanation goes here

D = params.D;
q = x2q(x);
tau = calctau(q, p, params);

f = D * tau;

end

