function [ yout , xout ] = sim_discrete( F , kfinal , x0 , y0 , u )
%sim_discrete: simulates a discrete system for kfinal timesteps, F is the
%discrete map function.
%   Detailed explanation goes here

yout = zeros( kfinal , length(y0) );
xout = zeros( kfinal , length(x0) );

yout(1,:) = y0';
xout(1,:) = x0';
for i = 2 : kfinal
    [ yplus , xplus ] = F( (i-1) , xout(i-1,:)' , u );
    yout(i,:) = yplus';
    xout(i,:) = xplus';
end

end


