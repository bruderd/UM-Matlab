% Makerange resmples data into a scale defined from 0 --> x where
% x is given as an input (range).
% 
% i.e. makerange(data, 100) is equivilant to makepercent because it
% resamples the data into tine units between 0 - 100
% 
% Ensure that first row is the time vector
% 
% Tim Dorn
% 12/12/2006
% ---------------------------------------------------


function new = makerange(old, range)

[tmp, l] = size(old);

t = old(1,:);
first = t(1);
last = t(end);

t = ((t - first)/(last-first))*range;

old(1,:) = t;
new = old;