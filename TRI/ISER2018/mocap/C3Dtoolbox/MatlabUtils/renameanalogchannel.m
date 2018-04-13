% Rename analog channel
% Tim Dorn
% 4th June 2009
% 
% --------------------------------------------------------------------
% Usage: renameanalogchannel(C3Dfile, origChan1, newChan1, origChan2, newChan2, ...)
% --------------------------------------------------------------------
% 
% Inputs:   C3Dfile: the C3D key structure from getEvents
%           origChanX: string of original channel name (case sensitive)
%           origChanX: string of new channel name      (case sensitive)
% 
% Outputs:  N/A
% 
% Notes
% -----
% 
% Original C3D file is overwritten.
% 
% e.g. renameanalogchannel('myfile.c3d', 'CH00', 'Fx1', 'CH01', 'Fy1');
%       renames CH00 to Fx1 and CH01 to Fy1
% 
% ---------------------------------------------------------------------


function renameanalogchannel(C3Dfile, varargin)

usage = 'Usage: renameanalogchannel(C3Dfile, origChan1, newChan1, origChan2, newChan2, ...)';


% Set up some initial parameters and do some initial checks
% ---------------------------------------------------------
larg = length(varargin);
    
if nargin < 2 || mod(larg, 2) == 1
    disp(usage)
    return
end

itf = c3dserver();
openc3d(itf, 0, C3Dfile);


% Collate tasks & Rename Channel
% ------------------------------
numChanges = larg/2;
for i = 1:numChanges
   origChannel = varargin{2*(i-1)+1};    % original analog channel
   newChannel = varargin{2*(i-1)+2};    % new analog channel
   found = 0;
   nIndex = itf.GetParameterIndex('ANALOG', 'LABELS');
   l = itf.GetParameterLength(nIndex);
   
   for j = 0:l-1
       if strcmp(itf.GetParameterValue(nIndex, j), origChannel)
           nSuccess = itf.SetParameterValue(nIndex, j, newChannel);
           fprintf('[%s]:  Analog Channel [%s] renamed to [%s]\n', ...
               C3Dfile, origChannel, newChannel)
           found = 1;
           break
       end
   end
   
   if ~found
       fprintf('[%s]:  Analog Channel [%s] not found\n', C3Dfile, origChannel);
   end
   
end

% Save new C3D file
% -----------------
nRet = itf.SaveFile('', -1);
closec3d(itf);

