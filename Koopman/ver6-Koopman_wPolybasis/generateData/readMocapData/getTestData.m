function [TR, PS] = getTestData(testname, params)
%getTestData: Read the test data from the pressure regulators and mocap
%system, line them up, and save a raw datafile with name 'testname' in the
%folder called 'rawdataFiles'

TR = struct;
PS = struct;

%%  Read in mocap data, save it in 'endeff' struct
disp('Please choose c3d file with your motion capture data in it...')
mocap = getC3Ddata(params);     % read in raw mocap data from C3D file
endeff = mocap2endeff(mocap);   % convert raw mocap data into endeffector coordinates

% Remove all data points that have NaNs in them
NaNindex = ~any(isnan(endeff.x),2);
PS.t = endeff.t(NaNindex,:);
PS.x = endeff.x(NaNindex,:);

%% Read in pressure data, save it in 'TR' struct (NOTE: some stuff here is specific to 3free system)
disp('Please choose the corresponding csv file with input pressure data...')
[pressure_file,pressure_path] = uigetfile;  % open dialog box for selecting file
pressureData = csvread([pressure_path, pressure_file]);    % matrix with columns: [time, input signals (V), measured pressures (Pa)]

TR.t = pressureData(:,1);
TR.pcontrol = pressureData(:,2:4);  % in (V). ONLY VALID BECAUSE VALVES 1-3 used!
TR.pmonitor = pressureData(:,10:12);       % in (Pa). ONLY VALID BECAUSE VALVES 1-3 used!

%% Partition the pressure data
% 
% % find all of the times when new control pressure signal sent (last entry is when pressure drops back to zero)
% TR.tsteps = [1, TR.t(1)];
% for i = 2:length(TR.t)
%    if any( TR.pcontrol(i,:) ~= TR.pcontrol(i-1,:) )
%        TR.tsteps = [ TR.tsteps; i, TR.t(i) ];
%    end
% end
% 
% % take average pressure at each step
% for i = 1:length(TR.tsteps) - 1
%     TR.psteps(i,:) = nanmean( TR.pmonitor( TR.tsteps(i,1)+5:TR.tsteps(i+1,1)-5, : ) );    % excludes some points that are in the transition region
% end
% % convert pressure in psteps into units of Pa
% TR.psteps = V2Pa(TR.psteps); 
% TR.pin = Vin2Pa(TR.pin, params);

%% Line up mocap and pressure data

% helpful plot for figuring out time offset
figure
hold on
plot(endeff.t, endeff.x(:,3), 'k')  % z wrt time
plot(TR.t, TR.pcontrol(:,3))    % input to 3rd actuator wrt t
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

% Line up pressure data with mocap data (i.e. mocap is the reference)
TR.t = TR.t - PS.t0;    % line the pressure data up with mocap data in time
TR.u = interp1(TR.t, TR.pcontrol, PS.t, 'linear', 0);    % interpolate to get all data points at same times, for points outside domains set equal to zero

% % This code lines mocap data with pressure data (i.e. pressure is the reference)
% PS.t = endeff.t + PS.t0;  % add the offset to sync the start times
% PS.x  = interp1(PS.t, endeff.x, TR.t);    % interpolate to get everything on the same time scale
% PS.tsteps = TR.tsteps;  % now they are synchronized in time


%% Remove Nans from data (already took care of this further up)

% % Remove all the NaNs from xsteps
% xsteps_nonan = PS.xsteps;
% psteps_nonan = TR.psteps;
% pin_nonan = TR.pin;
% 
% psteps_nonan(~any(~isnan(xsteps_nonan), 2),:)=[];
% pin_nonan(~any(~isnan(xsteps_nonan), 2),:)=[];
% xsteps_nonan(~any(~isnan(xsteps_nonan), 2),:)=[];
% 
% % store in structs
% TR.psteps_nonan = psteps_nonan;
% TR.pin_nonan = pin_nonan;
% PS.xsteps_nonan = xsteps_nonan;


%% Remove all points where pressure constraint was violated
%   This section currently not needed but left in because it shows how to selectively remove points 
%
% % just remove points when constraint was violated on FREEs 1 and 2
% TR.psteps_safe = TR.psteps_nonan( all(TR.psteps_nonan(:,1:2) < ones(size(TR.psteps_nonan(:,1))) * params.pmax(1:2), 2), : );
% PS.xsteps_safe = PS.xsteps_nonan( all(TR.psteps_nonan(:,1:2) < ones(size(TR.psteps_nonan(:,1))) * params.pmax(1:2), 2), : );
% TR.pin_safe = TR.pin_nonan( all(TR.psteps_nonan(:,1:2) < ones(size(TR.psteps_nonan(:,1))) * params.pmax(1:2), 2), : );

%% save data in a .mat file with fields t,x,u

t = PS.t;
x = PS.x;
u = TR.u;

% save the file using the system name that was set in params
[unique_fname, change_detect] = auto_rename(['rawdataFiles', filesep, testname, '.mat'], '0');
save(unique_fname, 't', 'x', 'u');


end
