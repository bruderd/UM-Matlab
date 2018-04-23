function error = debugQPconstraints( TR, PS, params)
%debugQPconstraints - Summary of this function goes here
%   Detailed explanation goes here

D = params.D;

for i = 1:length(TR.psteps)
    x = PS.xsteps(i,:)';
    p = TR.psteps(i,:)';
    
    q = x2q(x);
    Jq = calcJq(q, params);
    
    % calculate the elastomer force
    felast = calcFelast( x, params );
    
    Aeq = D*Jq';    % equality constraint makes boundary minima infeasable
    beq = -(felast);
    
    error(i,:) =  (Aeq * p - beq)';
end

end

