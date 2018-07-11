function [felast, params] = getFelast_2dof(TR, PS, params)
%calcFelast - performs system identification to determine the cumulative
%elastomer contribution.z

%% Calculate the average offset of x5 and x6
params.x5_offset = mean(PS.xsteps_nonan(:,5));
params.x6_offset = mean(PS.xsteps_nonan(:,6));


%% calculate y

for i = 1 : length(TR.psteps_nonan)
%     elast = calcf(PS.xsteps_safe(i,:)', TR.psteps_safe(i,:)', params); % no load force included for now
    elast = calcf(PS.xsteps_nonan(i,:)', TR.pin_nonan(i,:)', params); % no load force included for now
    elast = [0; 0; elast(3); elast(4); 0; 0];     % because all other directions are constrained, we only care about fitting the forces in the z and phi directions. And if the other directions are 0, there should be no elastomer force component.
    y(i,:) = -elast';
end

% % calculate y using all points that are not NaNs
% for i = 1 : length(psteps_nonan)
%     elast = calcf(xsteps_nonan(i,:)', psteps_nonan(i,:)', params); % no load force included for now
%     y(i,:) = -elast';
% end

%% fit polynomial

felast = struct;
% fit polynomial to elestomer force: felast
felast.x1 = MultiPolyRegress(PS.xsteps_nonan, y(:,1), 1);
felast.x2 = MultiPolyRegress(PS.xsteps_nonan, y(:,2), 1);
felast.x3 = MultiPolyRegress(PS.xsteps_nonan, y(:,3), 3, [0, 0, 3, 3, 0, 0], 'figure'); % I don't want elastomer force to be a function of x1, x2, x5, x6 since these are only nonzero due to sensor noise 
felast.x4 = MultiPolyRegress(PS.xsteps_nonan, y(:,4), 3, [0, 0, 3, 3, 0, 0], 'figure'); % I don't want elastomer force to be a function of x1, x2, x5, x6 since these are only nonzero due to sensor noise 
felast.x5 = MultiPolyRegress(PS.xsteps_nonan, y(:,5), 1);
felast.x6 = MultiPolyRegress(PS.xsteps_nonan, y(:,6), 1);

% NOTE: I chose 3rd order polynomial fits because it generated a predicted
% feasable region that most closely resembles the measured data.

% save the elastomer force within the params struct
params.felast = felast;


end
