function [TR, PS] = getTestData(pControlfile, pDatafile ,params)
%calcFelast - performs system identification to determine the cumulative
%elastomer contribution.z

%%  Read in mocap data, save it in 'endeff' struct
mocap = getC3Ddata(params);     % read in raw mocap data from C3D file
endeff = mocap2endeff(mocap);   % convert raw mocap data into endeffector coordinates

%% Read in pressure data, save it in 'TR' struct (NOTE: some stuff here is specific to 2dof system)
TR = struct;
TR.pin = csvread(pControlfile);    % the commanded voltages for each regulator
TR.pin = TR.pin(:,2:4); % remove first column of the data
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
    TR.psteps(i,:) = nanmean( TR.pmonitor( TR.tsteps(i,1)+5:TR.tsteps(i+1,1)-5, : ) );    % excludes some points that are in the transition region
end
% convert pressure in psteps into units of Pa
TR.psteps = V2Pa(TR.psteps); 
TR.pin = Vin2Pa(TR.pin, params);

%% Line up mocap and pressure data, and partition mocap data

% helpful plot for figuring out time offset
figure
hold on
plot(endeff.t, endeff.x(:,4), 'k')  % phi wrt time
stairs(TR.tsteps(1:end-1,2), TR.psteps(:,3)/10000)    % average pressure wrt time
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

%% Remove one outlier point (ONLY VALID FOR TEST 9!)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
TR.psteps = TR.psteps( find(PS.xsteps(:,3) > -0.007), : );
TR.pin = TR.pin( find(PS.xsteps(:,3) > -0.007), : );
PS.xsteps = PS.xsteps( find(PS.xsteps(:,3) > -0.007), : );


%% Remove Nans from data

% Remove all the NaNs from xsteps
xsteps_nonan = PS.xsteps;
psteps_nonan = TR.psteps;
pin_nonan = TR.pin;

psteps_nonan(~any(~isnan(xsteps_nonan), 2),:)=[];
pin_nonan(~any(~isnan(xsteps_nonan), 2),:)=[];
xsteps_nonan(~any(~isnan(xsteps_nonan), 2),:)=[];

% store in structs
TR.psteps_nonan = psteps_nonan;
TR.pin_nonan = pin_nonan;
PS.xsteps_nonan = xsteps_nonan;


%% Remove all points where pressure constraint was violated

% just remove points when constraint was violated on FREEs 1 and 2
TR.psteps_safe = TR.psteps_nonan( all(TR.psteps_nonan(:,1:2) < ones(size(TR.psteps_nonan(:,1))) * params.pmax(1:2), 2), : );
PS.xsteps_safe = PS.xsteps_nonan( all(TR.psteps_nonan(:,1:2) < ones(size(TR.psteps_nonan(:,1))) * params.pmax(1:2), 2), : );
TR.pin_safe = TR.pin_nonan( all(TR.psteps_nonan(:,1:2) < ones(size(TR.psteps_nonan(:,1))) * params.pmax(1:2), 2), : );


end
