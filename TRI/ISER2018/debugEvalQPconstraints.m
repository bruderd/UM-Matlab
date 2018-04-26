function const = debugEvalQPconstraints( TR, PS, params )
%evaluate the QP constraints
%   Detailed explanation goes here


for i = 1:length(PS.xsteps_nonan)
    x = PS.xsteps_nonan(i,:)';
    p = TR.psteps_nonan(i,:)';
    
    [H, f , A, b, Aeq, beq] = quadCost( x, params);
    
    model(i,:) = ( - A * [p; params.tol] + b )';    % constraint error. positive if constraint satisfied
    modelAll(i,:) = all( ( A * [p; params.tol] < b )' , 2);       % returns 0 for each constraint violated
    pressure(i,:) = (p < params.pmax');     % returns 0 if constraint violated

end

const.model = model;
const.modelAll = modelAll;  % 0 if any model constraint violated
const.pressure = pressure;