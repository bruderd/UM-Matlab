function tau = calctau( q, p, params )
%calcTau: Calculate the FREE forces, vertically concatenated into tau
%   Detailed explanation goes here

C = params.C;
Jq = calcJq(q, params);

tau = Jq'*p + C*q;

end

