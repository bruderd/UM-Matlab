% Makepercent resmples data into a percentage scale.
% It only uses the first row of the data, so the input
% can contain other rows of data, or simply a time vector.
% 
% Ensure that thie first row is in ascending order.
% 
% Tim Dorn
% 12/12/2006
% ---------------------------------------------------

function new = makepercent(old)

[tmp, l] = size(old);

t = old(1,:);
first = t(1);
last = t(end);

t = ((t - first)/(last-first))*100;

old(1,:) = t;
new = old;