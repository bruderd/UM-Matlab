% Extract raw EMG data EMG(t) from a C3D file
% 
% C3D --> EMG(t)
% 
% Tim Dorn
% June 2008
% 
% --------------------------------------------------------------------
% Usage: [EMGVecGlob, EMGset] = extractRawEMG(C3Dkey, emgSetName)
% --------------------------------------------------------------------
% 
% Inputs:   C3Dkey: the C3D key structure from getEvents
%           emgSetName: the label of EMG names contained in the EMG set
%               (this must be defined in loadlabels.m as glab.[emgSetName])
% 
% Outputs:  EMGVecGlob.MUSCNAME.time = time vector (starting from 0 in analog freq)
%           EMGVecGlob.MUSCNAME.data = raw EMG data vector (mV)
%           EMGVecGlob.MUSCNAME.name = string of muscle label
%           EMGVecGlob.MUSCNAME.event = events from C3Dkey for cropping
%           EMGVecGlob.MUSCNAME.c3dFile = C3D filename
%           EMGset = structure of cells containing information about the EMG labels
% 
% Notes
% -----
% 
% EMG LABELS ARE CONTAINED IN: loadLabels.m
% 
% ---------------------------------------------------------------------


function [EMGVecGlob, EMGset] = extractRawEMG(C3Dkey, emgSetName)

usage = 'Usage: [EMGVecGlob, EMGset] = extractRawEMG(C3Dkey, emgSetName)';


% Set up some initial parameters and do some initial checks
% ---------------------------------------------------------
if nargin ~= 2,
    disp(usage)
    return
end

if C3Dkey.allowed.EMG == 0
   error('Can not extract EMG due to a corrupted C3Dkey. Regenerate the C3Dkey and try again\n'); 
end


loadLabels;
EMGset = glab.(emgSetName);
numMuscles = length(EMGset);
mag = 1000;             % (from V to mV)


% Extract EMG data
% ----------------
itf = c3dserver();
openc3d(itf, 0, C3Dkey.c3dFile);
out.names = [];

for i = 1:numMuscles
    musc_name_orig = EMGset{i}{1};
    musc_name = removeSpaces(strtok(musc_name_orig, '.()[]{}'));
    out.names = [out.names, {musc_name}];
    data = mag*getanalogchannel(itf, musc_name_orig, ...
        C3Dkey.event.Vframes(1), C3Dkey.event.Vframes(end));
    time = C3Dkey.timeVec.Asec';
    
    % Slight reshape here to get correct number of frames
    % DONT NEED THIS SINCE I HAVE MODIFIED getanalogchannel.m (12/6/08)Tim
    %  data = data(:, 1:C3Dkey.numFrames.croppedA, :);

    EMGVecGlob.(musc_name).name = musc_name_orig;
    EMGVecGlob.(musc_name).data = double(data);
    EMGVecGlob.(musc_name).time = double(time);
    EMGVecGlob.(musc_name).event = C3Dkey.event;
    EMGVecGlob.(musc_name).c3dFile = C3Dkey.c3dFile;
end

closec3d(itf);

