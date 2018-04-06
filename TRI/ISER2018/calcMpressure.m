function [Pcontrol, feas, infeas] = calcMpressure(center, width, height, params)
%calcMpressure: Caluclates the input pressure for each point of the
%Michigan M
%   This is a quick way to set the dimensions of the M and caclulate the
%   input pressure at the same time. Alternatively you could generate
%   desired points and then run calcPressureTestPoints on those points

Mpoints = setMtestPoints( center, width, height );

feascount = 1;
infeascount = 1;
Pcontrol = zeros(length(Mpoints(:,1)), params.num);
for i = 1:length(Mpoints(:,1))
    [psol, exitflag] = calcPressure(Mpoints(i,:)', params);

    % show me which points are infeasable
    if exitflag > 0
        Pcontrol(i,:) = psol';
        feas(feascount,:) = Mpoints(i,3:4); % feasable
        feascount = feascount + 1;
    else
        infeas(infeascount,:) = Mpoints(i,3:4); % feasable
        infeascount = infeascount + 1;
    end
end

%% plot the feasable and infeasable points in different colors (optional)
figure
hold on
plot(feas(:,1), feas(:,2), 'b*')
plot(infeas(:,1), infeas(:,2), 'r*')
hold off

end