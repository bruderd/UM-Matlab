function [ f, gradf ] = costwgrad( x )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

f = cost(x);

gradf = gradcost(x);

end

