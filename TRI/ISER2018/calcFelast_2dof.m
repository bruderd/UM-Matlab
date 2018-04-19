function Felast = calcFelast_2dof(pDatafile ,params)
%calcFelast - performs system identification to determine the cumulative
%elastomer contribution.z

%%  Read in mocap data (comment out to save time for now CHANGE ME BACK DAN!)
mocap = getC3Ddata(params); 
endeff = mocap2endeff(mocap);

%% Read in pressure data
TR = struct;
pData = csvread(pDatafile);
TR.t = pData(:,1);
TR.pcontrol = pData(:,3:5);
TR.pmonitor = pData(:,11:13);

%% Do stuff with pressure data

% find all of the times when new control pressure signal sent (last entry is when pressure drops back to zero)
TR.tsteps = [1, TR.t(1)];
for i = 2:length(TR.t)
   if any( TR.pcontrol(i,:) ~= TR.pcontrol(i-1,:) )
       TR.tsteps = [ TR.tsteps; i, TR.t(i) ];
   end
end

% take average pressure at each step
for i = 1:length(TR.tsteps) - 1
    TR.psteps(i,:) = nanmean( TR.pmonitor( TR.tsteps(i,1)+5:TR.tsteps(i+1,1)-5, : ) );    % excludes some points that are in the transition region
end

%% Line up mocap and pressure data

% helpful plot for figuring out time offset
figure
hold on
plot(endeff.t, endeff.x(:,4), 'k')  % phi wrt time
stairs(TR.tsteps(1:end-1,2), TR.psteps)    % average pressure wrt time
hold off

t0_mocap = 0;   % the start of the mocap trial expressed in pressure time (figure this out via observation)
t_mocap = endeff.t + t0_mocap;  % add the offset to sync the start times
x_mocap  = interp1(t_mocap, endeff.x, TR.t);    % interpolate to get everything on the same time scale BUT REMEMBER YOU CAN'T INTERPOLATE EULER ANGLES SO YOU MUST CHANGE THIS!

% take average x at each step
for i = 1:length(TR.tsteps) - 1
    xsteps(i,:) = nanmean( x_mocap( TR.tsteps(i,1)+5:TR.tsteps(i+1,1)-5, : ) );    % excludes some points that are in the transition region
end

% plot to see if averaging x at each time step worked
figure
hold on
plot(TR.t, x_mocap(:,4), 'k')  % 
stairs(TR.tsteps(1:end-1,2), xsteps(:,4))    % 
hold off

%% Calculate elastomer force (y) at each point
for i = 1 : length(TR.tsteps)-1
    felast = calcf(xsteps(i,:)', TR.psteps(i,:)', params);
    y(i,:) = -felast';
end

% fit polynomial to elestomer force: e(x)
ex1 = MultiPolyRegress(xsteps, y(:,1), 2);
ex2 = MultiPolyRegress(xsteps, y(:,2), 2);
ex3 = MultiPolyRegress(xsteps, y(:,3), 2);
ex4 = MultiPolyRegress(xsteps, y(:,4), 2);
ex5 = MultiPolyRegress(xsteps, y(:,5), 2);
ex6 = MultiPolyRegress(xsteps, y(:,6), 2);


end
