% Create an event in a C3D file and save the C3D file
% 
% Tim Dorn
% June 2009
% 
% --------------------------------------------------------------------
% Usage: createEvent(c3dFile, foot, label, frame)
% --------------------------------------------------------------------
% 
% Inputs:   c3dFile: the name of the c3d file
% 
%           foot: 'R' for right foot event, 'L' for left foot event
%                 'G' for general foot event
% 
%           label: 'FS' for footstrike, 'FO' for footoff
%                  'GEN' for general event
% 
%           frame: video frame number to add the event at
% 
% 
% Outputs:  OVERWRITES the input c3d file with new events
%           (irreversible with this code so make sure to 
%           backup the original c3d file first!)
% 
% ---------------------------------------------------------------------

function createEvent(c3dFile, foot, label, frame)

usage = 'createEvent(c3dFile, foot, label, frame)';

if nargin ~= 4
    error(usage)
end

switch upper(foot)
    case 'R'
        contexts = 'Right';
    case 'L'
        contexts = 'Left';
    case 'G'
        contexts = 'General';         
    otherwise
        error(usage)
end

switch upper(label)
    case 'FS'
        icon_ids = 1;
        labels = 'Foot Strike';
        descriptions = labels;
        generic_flags = 0;
    case 'FO'
        icon_ids = 2;
        labels = 'Foot Off';
        descriptions = labels;
        generic_flags = 0;
    case 'GEN'
        icon_ids = 0;
        labels = 'General Event';
        descriptions = labels;
        generic_flags = 1;        
    otherwise
        error(usage)
end


% ADD EVENT
% ----------
itf = c3dserver();
openc3d(itf, 0, c3dFile);


% check if events group exists. if not, then create the group
% and all parameters that come under it.
eventGroupName = 'EVENT';      % do not change
eventIndex = itf.GetGroupIndex(eventGroupName);
if eventIndex < 0
    itf.AddGroup(0, eventGroupName, 'Event Labels', char(0));
    itf.AddParameter('USED', 'Used', eventGroupName,                   char(0), 2, 0, int16([1,0]), {});
    itf.AddParameter('CONTEXTS', 'Contexts', eventGroupName,           char(0), -1, 2, int16([16,0]), {});
    itf.AddParameter('ICON_IDS', 'Icon IDs', eventGroupName,           char(0), 2, 1, int16([1,0]), {});
    itf.AddParameter('LABELS', 'Labels', eventGroupName,               char(0), -1, 2, int16([32,0]), {});
    itf.AddParameter('DESCRIPTIONS', 'Descriptions', eventGroupName,   char(0), -1, 2, int16([80,0]), {});
    itf.AddParameter('SUBJECTS', 'Subjects', eventGroupName,           char(0), -1, 2, int16([32,0]), {});
    itf.AddParameter('TIMES', 'Times', eventGroupName,                 char(0), 4, 2, int16([2,0]), {});
    itf.AddParameter('GENERIC_FLAGS', 'Generic Flags', eventGroupName, char(0), 1, 1, int16([1,0]), {});
    
    index3 = itf.GetParameterIndex(eventGroupName, 'ICON_IDS');
    itf.RemoveParameterData(index3, 0);
    index8 = itf.GetParameterIndex(eventGroupName, 'GENERIC_FLAGS');
    itf.RemoveParameterData(index8, 0);
end


index1 = itf.GetParameterIndex(eventGroupName, 'USED');
oldnumEvents = itf.GetParameterValue(index1, 0);
itf.SetParameterValue(index1, 0, oldnumEvents+1);

index2 = itf.GetParameterIndex(eventGroupName, 'CONTEXTS');
itf.AddParameterData(index2, 1);
itf.SetParameterValue(index2, oldnumEvents, contexts);
itf.GetParameterValue(index2, oldnumEvents);

index3 = itf.GetParameterIndex(eventGroupName, 'ICON_IDS');
itf.AddParameterData(index3, 1);
itf.SetParameterValue(index3, oldnumEvents, icon_ids);

index4 = itf.GetParameterIndex(eventGroupName, 'LABELS');
itf.AddParameterData(index4, 1);
itf.SetParameterValue(index4, oldnumEvents, labels);

index5 = itf.GetParameterIndex(eventGroupName, 'DESCRIPTIONS');
itf.AddParameterData(index5, 1);
itf.SetParameterValue(index5, oldnumEvents, descriptions);

index6 = itf.GetParameterIndex(eventGroupName, 'SUBJECTS');
if oldnumEvents == 0
    subjects = 'SUBJECT';
else
    subjects = itf.GetParameterValue(index6, oldnumEvents-1);
end
itf.AddParameterData(index6, 1);
itf.SetParameterValue(index6, oldnumEvents, subjects);

index7 = itf.GetParameterIndex(eventGroupName, 'TIMES');
times = (frame-1)/itf.GetVideoFrameRate;
itf.AddParameterData(index7, 1);
itf.SetParameterValue(index7, 2*oldnumEvents+1, times);

index8 = itf.GetParameterIndex(eventGroupName, 'GENERIC_FLAGS');
itf.AddParameterData(index8, 1);
itf.SetParameterValue(index8, oldnumEvents, generic_flags);




% SAVE C3D FILE
% --------------

nRet = itf.SaveFile('', -1);
if nRet ~= 1
    fprintf('WARNING: C3D file could not be saved!\n');
end

closec3d(itf);
