% test the calcPressure function. Outputs a plot of the exit flags over a
% range of x values for that you can see where the system is infeasable and
% such. This is written to by used with the 2 DOF rotation and translation rig
clear;

params = setParams;

dl = -0.03:0.005:0.03;
dphi = -(3*pi):0.2:(3*pi);

testPoints = combvec(dl,dphi);

% exitflag = zeros(1, length(testPoints));
% for i = 1:length(testPoints)
%    dli = testPoints(1,i);
%    dphii = testPoints(2,i);
%    [~, exitflag(i)] = calcPressure([0,0,dli,0,0,dphii]', params); 
% end
% [X,Y] = meshgrid(testPoints(1,:), testPoints(2,:));

[X,Y] = meshgrid(dl, dphi);

Z = zeros(size(X));
for i = 1:length(X(:,1))
    for j = 1:length(X(1,:))
        dli = X(i,j);
        dphii = Y(i,j);
        [~, Z(i,j)] = calcPressure([0,0,dli,0,0,dphii]', params);
    end
end

%% plot it
figure
pcolor(X,Y,Z);
colorbar
