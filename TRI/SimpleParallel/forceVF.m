% forceVF.m
% Plot the force vector [F, M]' over range of pressures

% set the pressure step sizes based on the min/max from setParams.m
dP_l = (params.Prange_l(2) - params.Prange_l(1)) / params.steps_l;
dP_r = (params.Prange_r(2) - params.Prange_r(1)) / params.steps_r;

% initializations
P_l = zeros(1,params.steps_l);
P_r = zeros(1,params.steps_r);
F = zeros(params.steps_l, params.steps_r);
M = zeros(params.steps_l, params.steps_r);


for i = 1:params.steps_l
    
    P_l(i) = i * dP_l;
    
    for j = 1:params.steps_r
        
        P_r(j) = j * dP_r;
             
        netForce = netF([0;1], [P_l(i); P_r(j)], params);
        F(i,j) = netForce(1);
        M(i,j) = netForce(2);
        
    end
    
end

P_kPa = [P_l; P_r] * 1e-3;
[X,Y] = meshgrid(P_kPa(1,:), P_kPa(2,:));

figure
forcefield = quiver(Y,X,F,M*10);