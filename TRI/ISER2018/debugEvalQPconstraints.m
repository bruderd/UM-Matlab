function const = debugEvalQPconstraints( TR, PS, params )
%evaluate the QP constraints
%   Detailed explanation goes here


for i = 1:length(PS.xsteps_nonan)
    const.x(i,:) = PS.xsteps_nonan(i,:);
    x = const.x(i,:)';
    const.p(i,:) = TR.psteps_nonan(i,:)';
    p = const.p(i,:)';
    
    [H, f , A, b, Aeq, beq] = quadCost( x, params);
    
    model(i,:) = ( - A * [p; params.tol'] + b )';    % constraint error. positive if constraint satisfied
    modelFeas(i,:) = ( A * [p; params.tol'] < b )';
    modelAll(i,:) = all( ( A * [p; params.tol'] < b )' , 2);       % returns 0 for each constraint violated
    pressure(i,:) = (p < params.pmax');     % returns 0 if constraint violated

end

const.model = model;
const.modelFeas = modelFeas;
const.modelAll = modelAll;  % 0 if any model constraint violated
const.pressure = pressure;