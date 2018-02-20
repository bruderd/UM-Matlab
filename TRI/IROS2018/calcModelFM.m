function predForces = calcModelFM( testPoints, params )
%calcModelFM: Predict the force and moment at each test point base on model
%   Detailed explanation goes here

predForces = zeros(size(testPoints,1), 2);
for i = 1:size(testPoints,1)
    P = testPoints(i, 3:end)';
    x = [0,0, testPoints(i,1)*10^(-3), 0, 0, deg2rad(testPoints(i,2))]';
    zeta = calczeta(x, P, params);
    predForces(i,:) = [zeta(3), zeta(6)];
end

end

