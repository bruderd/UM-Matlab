function [felast, params] = getFelast_2dof(TR, PS, params)
%calcFelast - performs system identification to determine the cumulative
%elastomer contribution.z

% Remove all the NaNs from xsteps
xsteps_nonan = PS.xsteps;
psteps_nonan = TR.psteps;
psteps_nonan(~any(~isnan(xsteps_nonan), 2),:)=[];
xsteps_nonan(~any(~isnan(xsteps_nonan), 2),:)=[];

for i = 1 : length(psteps_nonan)
    elast = calcf(xsteps_nonan(i,:)', psteps_nonan(i,:)', params); % no load force included for now
    y(i,:) = -elast';
end

felast = struct;
% fit polynomial to elestomer force: felast
felast.x1 = MultiPolyRegress(xsteps_nonan, y(:,1), 2);
felast.x2 = MultiPolyRegress(xsteps_nonan, y(:,2), 2);
felast.x3 = MultiPolyRegress(xsteps_nonan, y(:,3), 2);
felast.x4 = MultiPolyRegress(xsteps_nonan, y(:,4), 2);
felast.x5 = MultiPolyRegress(xsteps_nonan, y(:,5), 2);
felast.x6 = MultiPolyRegress(xsteps_nonan, y(:,6), 2);

% save the elastomer force within the params struct
params.felast = felast;

end
