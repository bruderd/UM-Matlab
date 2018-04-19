function Felast = calcFelast_2dof(pDatafile ,params)
%calcFelast - performs system identification to determine the cumulative
%elastomer contribution.z

%%  Read in mocap data (comment out to save time for now CHANGE ME BACK DAN!)
% mocap = getC3Ddata(params); 
% endeff = mocap2endeff(mocap);

%% Read in pressure data
TR = struct;
pData = csvread(pDatafile);
TR.t = pData(:,1);
TR.pcontrol = pData(:,3:5);
TR.pmonitor = pData(:,11:13);

%% Do stuff with pressure data

% find all of the times when new control pressure signal sent
TR.tsteps = [1, TR.t(1)];
for i = 2:length(TR.t)
   if any( TR.pcontrol(i,:) ~= TR.pcontrol(i-1,:) )
       TR.tsteps = [ TR.tsteps; i, TR.t(i) ];
   end
end


% take average pressure at each pressure configuration
for i = 1:length(TR.tsteps) - 1
    ppoint(i,:) = mean( TR.pmonitor( TR.tsteps(i,1)+5:TR.tsteps(i+1,1)-5, : ) );    % excludes some points that are in the transition region
end

end
