function [y, ud, yd] = test_narx( valdata )
% test out the neural net model

% preallocation
TS = length(valdata.t);
y = zeros(TS,2);
ud = zeros(TS,3);
yd = zeros(TS,2);

% initial values
ud(1:2,:) = valdata.u(1:2,:);
yd(1:2,:) = valdata.x(1:2,:);
y(1:2,:) = valdata.x(1:2,:);

for i = 2 : TS-1
    
    [y1,xf1,xf2] = larm_narx_1del( ud(i,:)' , yd(i,:)' , ud(i-1,:)' , yd(i-1,:) );
    
    y(i+1,:) = y1';
    ud(i+1,:) = xf1';
    yd(i+1,:) = xf2';
end