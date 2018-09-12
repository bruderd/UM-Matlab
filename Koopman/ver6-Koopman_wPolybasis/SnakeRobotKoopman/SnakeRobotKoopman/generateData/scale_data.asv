function [ x_scaled, u_scaled ] = scale_data( x , u )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

n = size(x,2);
p = size(u,2);

% maximum value for each of the states (doesn't include input), in a row vector
maxStates = max( abs( max(x, [],1) ) , abs( min(x, [],1) ));

% scale each state so that it never goes outside the range [-0.1, 0.1]
x_scaled = x ./ (10 * maxStates);
% x_scaled = x;

% scale the input so that it never goes above 0.1
u_scaled = u ./ 100;

end

