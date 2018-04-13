% Generate a marker *.trc file readable by OpenSim
% 
% Tim Dorn
% June 2009
% 
% --------------------------------------------------------------------
% Usage: generateTrcFile(C3Dkey, markerpos, markerset)
% --------------------------------------------------------------------
% 
% Inputs:   C3Dkey: the C3D key structure from getEvents
%           markerpos = array of marker positions
%                   for M markers: should contain 1+3M columns
%                   (time + XYZ of each marker)
%           markerset = cell array of strings containing the names of markers
%                   e.g. markerset = {'M1', 'M2', 'M3'};
% 
% Outputs:  output trc file
% 
% 
% Important note: the opensim convention is to output the GRF and CoP from
% start time to end time (columns) and RIGHT foot forces/positions then
% LEFT foot forces/positions (rows).
% 
% ----------------------------------------------------------------------

function generateTrcFile(C3Dkey, markerpos, markerset)
    
PathFileType = 4;
name = sprintf('%s.trc', C3Dkey.c3dFile);
datatype = '(X/Y/Z)';
DataRate = C3Dkey.vFreq;
CameraRate = C3Dkey.vFreq;
NumFrames = length(markerpos(:,1));
NumMarkers = length(markerset);
Units = 'mm';
OrigDataRate = C3Dkey.vFreq;
OrigDataStartFrame = C3Dkey.event.Vframes(1);
OrigNumFrames = NumFrames;

% Adjust frames & times so that the effective trial period
% starts at frame 1 at time 0
markerpos(:,1) = markerpos(:,1) - C3Dkey.event.Vframes(1) + 1;
t = C3Dkey.timeVec.Vsec';
frame = markerpos(:,1);
markerpos = [t, markerpos(:,2:end)];



% TRC File Header
% ---------------

fid = fopen(name, 'w');
if fid < 0
    fprintf('\nERROR: %s could not be opened for writing...\n\n', name);
    return
end
fprintf(fid, 'PathFileType\t%d\t%s\t%s\t\n', PathFileType, datatype, name);
fprintf(fid, 'DataRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\tOrigDataRate\tOrigDataStartFrame\tOrigNumFrames\n');
fprintf(fid, '%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\n', ...
    DataRate, CameraRate, NumFrames, NumMarkers, Units, OrigDataRate, OrigDataStartFrame, OrigNumFrames);
fprintf(fid, 'Frame#\tTime\t');


% TRC File Body
% -------------

for i = 1:NumMarkers
    fprintf(fid, '%s\t\t\t', markerset{i});
end
fprintf(fid, '\n\t\t');

for i = 1:NumMarkers
    fprintf(fid, 'X%d\tY%d\tZ%d\t', i, i, i);
end
fprintf(fid, '\n\n');

% marker position values
for i = 1:NumFrames
    fprintf(fid, '%d\t', frame(i));
    fprintf(fid, '%.5f\t', markerpos(i,:));
    fprintf(fid, '\n');
end

fclose(fid);
fprintf('Saved (tab delimited) marker positions to: %s\n', name);

