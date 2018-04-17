function [ endeff ] = mocap2endeff( mocap )
%mocap2x: Uses mocap marker position data to compute the state of the end
%   Detailed explanation goes here

% calculate center of top block as the average position of top block markers over whole trial
total = zeros(1,3);
for i = 1 : size(mocap.topData, 3)
   total = total + nanmean(mocap.topData(:,:,i));
   topcent = total / i;
end





% create output struct with state and timestamp of end effector
endeff = struct;
endeff.x = x;
endeff.t = t;

end

