function [ mocap, frames ] = getC3Ddata(params, C3Dfile)
%getC3Ddata - extracts marker positions from a c3d data file
%   Detailed explanation goes here

topIDs = params.topIDs;
effIDs = params.effIDs;

itf = c3dserver();

if exist('C3Dfile','var')
    openc3d(itf, varargin); % opens c3d file whose path is specified by optional second argument
else
    openc3d(itf, 1);    % opens a dialog box where you can choose c3d file
end

frames = nframes(itf);
time = (1/60) * [1:frames]';

% get the data from the C3D file
xyz = get3dtargets(itf, 0); % second argument determines if residual should be included

topData = zeros(frames, 3, length(topIDs));
for i = 1 : length(topIDs)
    topLEDs(i) = string( strcat( 'xyz.M', sprintf('%03d', topIDs(i)) ) );
    topData(:,:,i) = eval( strcat( 'xyz.M', sprintf('%03d', topIDs(i)) ) );
end

effData = zeros(frames, 3, length(effIDs));
for i = 1 : length(effIDs)
    effLEDs(i) = string( strcat( 'xyz.M', sprintf('%03d', effIDs(i)) ) );
    effData(:,:,i) = eval( strcat( 'xyz.M', sprintf('%03d', effIDs(i)) ) );
end

mocap = struct;
mocap.xyz = xyz;    % raw mocap position data for each marker
mocap.t = time;     % time of each sample (s)
mocap.topLEDs = topLEDs;    % string array with the names of the LED markers for the top block
mocap.endeffLEDs = effLEDs; % string array with the names of the LED markers for the end effector
mocap.topData = topData;    % 3D matrix containing the marker data for all top block markers
mocap.effData = effData;    % 3D matrix containing the marker data for all end effector markers

end

