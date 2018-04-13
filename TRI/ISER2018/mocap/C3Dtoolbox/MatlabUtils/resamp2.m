% Resample Data
% Tim Dorn
% 20/7/2006
% 
% Usage: sampled = resamp2(raw, frames)
% 
% Use this to resample a set of data into a certain number
% of frames... Unlike resamp1, this does not require a time
% vector at the start
% ---------------------------------------------------------

function sampled = resamp2(raw, frames)

[m,n] = size(raw);

tr = 0:n-1;
ts = makerange(0:frames-1, n-1);

sampled = [];
for i = 1:m
    sampled = [sampled ; interp1(tr, raw(i,:), ts, 'spline')];
end

ts = 0:frames-1;

% close all
% for i = 1:m
%     figure
%     plot(tr, raw(i,:), 'b')
%     hold on
%     plot(ts, sampled(i,:), 'r')
% end
