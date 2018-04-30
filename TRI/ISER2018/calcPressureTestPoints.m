function [Pcontrol, feas, infeas] = calcPressureTestPoints(testPoints, params)
%calcMpressure: Caluclates the input pressure for each testPoint
%   Width of testPoints is 6 since rows are the 6 dimensional test points,
%   height is the total number of test points.

feascount = 1;
infeascount = 1;
Pcontrol = zeros(length(testPoints(:,1)) + 1, params.num);

Pcontrol(1,:) = zeros(1, params.num);   % for calibration, first point is zero pressure

for i = 1:length(testPoints(:,1))
    [psol, exitflag] = calcPressure(testPoints(i,:)', params);
    
    % show me which points are infeasable
    if exitflag > 0
        Pcontrol(i+1,:) = psol(1:params.num)';    % just get the pressure values not the tolerance
        feas(feascount,:) = testPoints(i,3:4); % feasable
        feascount = feascount + 1;
    else
        Pcontrol(i+1,:) = zeros(1, params.num);    % set pressure = 0 at infeasable points
        infeas(infeascount,:) = testPoints(i,3:4); % infeasable
        infeascount = infeascount + 1;
    end
end

%% plot the feasable and infeasable points in different colors (optional)
figure
hold on
if exist('feas', 'var')
    plot(feas(:,1), feas(:,2), 'b*');
end
if exist('infeas', 'var')
    plot(infeas(:,1), infeas(:,2), 'r*')
end
hold off

end