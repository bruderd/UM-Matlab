function Felast = calcFelast(pDatafile ,params)
%calcFelast - performs system identification to determine the cumulative
%elastomer contribution.z


mocap = getC3Ddata(params); 
endeff = mocap2endeff(mocap);

% find all of the times when new control pressure signal sent
TR = struct;
pData = csvread('pDatafile');
TR.t = pData(:,1);
TR.pcontrol = pData(:,2:9);
TR.pmonitor = pData(:,10:17);
