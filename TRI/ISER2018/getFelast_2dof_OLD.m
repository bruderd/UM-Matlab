function [felast, params] = getFelast_2dof_OLD(TR, PS, params)

%%  Read in mocap data, save it in 'endeff' struct
mocap = getC3Ddata(params);     % read in raw mocap data from C3D file
endeff = mocap2endeff(mocap);   % convert raw mocap data into endeffector coordinates

%% Read in pressure data, save it in 'TR' struct
TR = struct;
pData = csvread(pDatafile);
TR.t = pData(:,1);
TR.pcontrol = pData(:,3:5);
TR.pmonitor = pData(:,11:13);

%% Partition the pressure data

% find all of the times when new control pressure signal sent (last entry is when pressure drops back to zero)
TR.tsteps = [1, TR.t(1)];
for i = 2:length(TR.t)
   if any( TR.pcontrol(i,:) ~= TR.pcontrol(i-1,:) )
       TR.tsteps = [ TR.tsteps; i, TR.t(i) ];
   end
end

% take average pressure at each step
for i = 1:length(TR.tsteps) - 1
    TR.psteps_V(i,:) = nanmean( TR.pmonitor( TR.tsteps(i,1)+5:TR.tsteps(i+1,1)-5, : ) );    % excludes some points that are in the transition region
    TR.pcontrolsteps_V(i,:) = nanmean( TR.pcontrol( TR.tsteps(i,1)+5:TR.tsteps(i+1,1)-5, : ) );    % excludes some points that are in the transition region
end
TR.psteps = TR.psteps_V * (150/10) * (6894.76/1);    % convert pressure voltage back into psi
TR.pcontrolsteps = TR.pcontrolsteps_V * (params.TRpsimax/10) * (6894.76/1); 

%% Line up mocap and pressure data, and partition mocap data

% helpful plot for figuring out time offset
figure
hold on
plot(endeff.t, endeff.x(:,4), 'k')  % phi wrt time
stairs(TR.tsteps(1:end-1,2), TR.psteps_V(:,3))    % average pressure wrt time
title('Mocap and Pressure Data Together')
legend('Mocap: PhaseSpace', 'Pressure: TR')
hold off

% Prompt the user to identify the time offset between the data sets from
% the plot
prompt = 'What is the time offset between the PhaseSpace and TR data, in seconds?:';
offset = input(prompt);     % the start of the mocap trial expressed in pressure time
if isnumeric(offset)
    PS.t0 = offset;
else
    disp('Invalid. You must enter a numeric value.');
    return;
end

PS.t = endeff.t + PS.t0;  % add the offset to sync the start times
PS.x  = interp1(PS.t, endeff.x, TR.t);    % interpolate to get everything on the same time scale BUT REMEMBER YOU CAN'T INTERPOLATE EULER ANGLES SO YOU MUST CHANGE THIS!
PS.tsteps = TR.tsteps;  % now they are synchronized in time

% take average x at each step
for i = 1:length(TR.tsteps) - 1
    PS.xsteps(i,:) = nanmean( PS.x( TR.tsteps(i,1)+5:TR.tsteps(i+1,1)-5, : ) );    % excludes some points that are in the transition region
end

% plot to see if averaging x at each time step worked (for DEBUGGING)
figure
hold on
plot(TR.t, PS.x(:,4), 'k')  % 
stairs(TR.tsteps(1:end-1,2), PS.xsteps(:,4))    % 
hold off

%% Elastomer polynomial fit

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
