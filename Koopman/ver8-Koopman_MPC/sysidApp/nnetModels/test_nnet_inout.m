function [ zeta, y ] = test_nnet_inout( valdata )
% test out the neural net model

% preallocation and assigment
TS = length(valdata.t);
zeta = zeros(TS,5);
y = zeros(TS,2);
u = [ valdata.u ];   % add one extra row
yreal = valdata.x;

% initial values
y(1,:) = valdata.x(2,:);
zeta(1:2,:) = [ yreal(1:2,:) , u(1:2,:) ];

for i = 2 : TS
    
    zeta(i,:) = [ y(i-1,:) , u(i,:) ];

    [yout,xf1] = larm_nnet_inout_1del( zeta(i,:)' , zeta(i-1,:)' );
    
    y(i,:) = yout';
%     zeta(i+1,:) = [ y(i,:) , u(i+1,:), u(i+2,:) ];

end