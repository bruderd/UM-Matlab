function [felast, params] = getFelast_2dof(TR, PS, params)
%calcFelast - performs system identification to determine the cumulative
%elastomer contribution.z

% Remove all the NaNs from xsteps
xsteps_nonan = PS.xsteps;
psteps_nonan = TR.psteps;
psteps_nonan(~any(~isnan(xsteps_nonan), 2),:)=[];
xsteps_nonan(~any(~isnan(xsteps_nonan), 2),:)=[];

% calculate y using only points where pressure constraint not violated
for i = 1 : length(TR.psteps_safe)
    elast = calcf(PS.xsteps_safe(i,:)', TR.psteps_safe(i,:)', params); % no load force included for now
    y(i,:) = -elast';
end

% % calculate y using all points that are not NaNs
% for i = 1 : length(psteps_nonan)
%     elast = calcf(xsteps_nonan(i,:)', psteps_nonan(i,:)', params); % no load force included for now
%     y(i,:) = -elast';
% end

felast = struct;
% fit polynomial to elestomer force: felast
felast.x1 = MultiPolyRegress(PS.xsteps_safe, y(:,1), 1);
felast.x2 = MultiPolyRegress(PS.xsteps_safe, y(:,2), 1);
felast.x3 = MultiPolyRegress(PS.xsteps_safe, y(:,3), 3);
felast.x4 = MultiPolyRegress(PS.xsteps_safe, y(:,4), 3);
felast.x5 = MultiPolyRegress(PS.xsteps_safe, y(:,5), 3);
felast.x6 = MultiPolyRegress(PS.xsteps_safe, y(:,6), 3);

% save the elastomer force within the params struct
params.felast = felast;

end
