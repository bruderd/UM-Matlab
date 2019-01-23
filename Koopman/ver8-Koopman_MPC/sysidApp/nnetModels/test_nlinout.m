function [y, ud] = test_nlinout( valdata )
% test out the neural net model

% preallocation
TS = length(valdata.t);
y = zeros(TS,2);
ud = valdata.u; %zeros(TS,3);

% initial values
ud(1:2,:) = valdata.u(1:2,:);
y(1:2,:) = valdata.x(1:2,:);

for i = 2 : TS-1
    
    [y1,xf1] = larm_nlinout_1del( ud(i,:)' , ud(i-1,:)' );
    
    y(i+1,:) = y1';
%     ud(i+1,:) = xf1';
end