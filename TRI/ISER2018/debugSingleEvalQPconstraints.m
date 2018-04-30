function const = debugSingleEvalQPconstraints( p, x, params )
%evaluate the QP constraints at a single point
%   Detailed explanation goes here


for i = 1:1 % loop unneccesary. Residue of a copy/past hack job
    
    [H, f , A, b, Aeq, beq] = quadCost( x, params);
    
    Amatrix(:,:,i) = A;
    model(i,:) = ( - A * [p; params.tol'] + b )';    % constraint error. positive if constraint satisfied
    modelFeas(i,:) = ( A * [p; params.tol'] < b )';
    modelAll(i,:) = all( ( A * [p; params.tol'] < b )' , 2);       % returns 0 for each constraint violated
    pressure(i,:) = (p <= params.pmax');     % returns 0 if constraint violated
    pressure2(i,:) = (p >= params.pmin');     % returns 0 if constraint violated

end

const.model = model;
const.modelFeas = modelFeas;
const.modelAll = modelAll;  % 0 if any model constraint violated
const.pressure = pressure;
const.pressure2 = pressure2;
const.Amatrix = Amatrix;