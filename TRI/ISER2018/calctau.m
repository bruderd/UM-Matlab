function tau = calctau( q, p, params )
%calcTau: Calculate the FREE fiber forces, vertically concatenated into tau
%   Detailed explanation goes here

C = params.C;
Jq = calcJq(q, params);

tau = Jq'*p;

% tau = Jq'*p + C*q;    % with linear elastomer force

end

