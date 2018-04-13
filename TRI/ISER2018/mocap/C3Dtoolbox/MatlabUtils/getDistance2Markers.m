% Get distance between two markers
% Tim Dorn
% April 2009
% 
% --------------------------------------------------------------------
% Usage: [distXYZ, distNORM, distNORMAVR] = getDistance2Markers(markerStruct, M1, M2)
% --------------------------------------------------------------------
% 
% Inputs:   markerStruct: the marker structure from getMarkers.m
%           M1: string of 1st marker name 
%           M1: string of 1st marker name 
% 
% 
% Outputs:  distXYZ: XYZ distances from M1 to M2 for each time frame
%           distAVR: average of distXYZ over all time frames
%           distAVRNORM: average norm of distXYZ over all time frames
% 
% 
% Notes:    N/A
% 
% --------------------------------------------------------------------

function [distXYZ, distAVR, distAVRNORM, midpt] = getDistance2Markers(markerStruct, M1, M2)


% find marker position in structure
M1_1 = findMarkerPos(M1, markerStruct.label);
M1_pos = [M1_1, M1_1+1, M1_1+2];

if nargin == 2,
    distXYZ = markerStruct.data(:, M1_pos);
elseif nargin == 3,
    M2_1 = findMarkerPos(M2, markerStruct.label);
    M2_pos = [M2_1, M2_1+1, M2_1+2];
    distXYZ = markerStruct.data(:, M2_pos) - markerStruct.data(:, M1_pos);
    midpt = (mean(markerStruct.data(:, M2_pos)) + mean(markerStruct.data(:, M1_pos))) / 2;
end

distAVR = mean(distXYZ, 1);
distAVRNORM = norm(distAVR);




% Subfunction to find marker string in the
% list of labels and return the column index
% ------------------------------------------
function colIndex = findMarkerPos(M, labels)

l = length(labels);
for i = 1:l
    if strcmp(M, strtok(labels(i), '-'))
        colIndex = i;
        return;
    end
    
end
colIndex = 0;
fprintf('No matching marker found in label list...\n')
return;


