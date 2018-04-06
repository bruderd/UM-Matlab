function [Pcontrol, feas, infeas] = calcPressureTestPoints(testPoints, params)
%calcMpressure: Caluclates the input pressure for each testPoint
%   Width of testPoints is 6 since rows are the 6 dimensional test points,
%   height is the total number of test points.

feascount = 1;
infeascount = 1;
Pcontrol = zeros(length(testPoints(:,1)), params.num);
for i = 1:length(testPoints(:,1))
    [psol, exitflag] = calcPressure(testPoints(i,:)', params);
    Pcontrol(i,:) = psol';
    
    % show me which points are infeasable
    if exitflag > 0
        feas(feascount,:) = testPoints(i,3:4); % feasable
        feascount = feascount + 1;
    else
        infeas(infeascount,:) = testPoints(i,3:4); % feasable
        infeascount = infeascount + 1;
    end
end

%% plot the feasable and infeasable points in different colors (optional)
% figure
% hold on
% plot(feas(:,1), feas(:,2), 'b*')
% plot(infeas(:,1), infeas(:,2), 'r*')
% hold off

end