function [felast, params] = getFelast_2dof(TR, PS, params)
%calcFelast - performs system identification to determine the cumulative
%elastomer contribution.z


% calculate y using only points where pressure constraint not violated
for i = 1 : length(TR.psteps_nonan)
%     elast = calcf(PS.xsteps_safe(i,:)', TR.psteps_safe(i,:)', params); % no load force included for now
    elast = calcf(PS.xsteps_nonan(i,:)', TR.pin_nonan(i,:)', params); % no load force included for now
    elast = [0; 0; elast(3); elast(4); 0; 0];     % because all other directions are constrained, we only care about fitting the forces in the z and phi directions
    y(i,:) = -elast';
end

% % calculate y using all points that are not NaNs
% for i = 1 : length(psteps_nonan)
%     elast = calcf(xsteps_nonan(i,:)', psteps_nonan(i,:)', params); % no load force included for now
%     y(i,:) = -elast';
% end

felast = struct;
% fit polynomial to elestomer force: felast
felast.x1 = MultiPolyRegress(PS.xsteps_nonan, y(:,1), 1);
felast.x2 = MultiPolyRegress(PS.xsteps_nonan, y(:,2), 1);
felast.x3 = MultiPolyRegress(PS.xsteps_nonan, y(:,3), 3, 'figure');
felast.x4 = MultiPolyRegress(PS.xsteps_nonan, y(:,4), 3, 'figure');
felast.x5 = MultiPolyRegress(PS.xsteps_nonan, y(:,5), 1, 'figure');
felast.x6 = MultiPolyRegress(PS.xsteps_nonan, y(:,6), 1, 'figure');

% save the elastomer force within the params struct
params.felast = felast;

end
