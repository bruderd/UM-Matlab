% Batch Process multiple stride EMG signals from a single C3D file
% Tim Dorn
% July 2009
% 
% --------------------------------------------------------------------
% Usage: [eVecGlob, EMGVecGlob] = multipleEMGprocess(c3dFile, emgSetName, emgProcessTasks, interval4Time)
% --------------------------------------------------------------------
% 
% Inputs:   c3dFile = the name of the c3d file
% 
%           emgSetName: the label of EMG names contained in the EMG set
%               (this must be defined in loadlabels.m as a cell: glab.[emgSetName])
% 
%           emgProcessTasks: the EMG processing options in the order of execution
%               (this must be defined in loadlabels.m as a cell: glab.[emgProcessTasks])
% 
%           interval4Time = the interval number to set for the time vector
% 
% 
% Outputs:  eVecGlob = structure of processed EMG
%           EMGVecGlob = structure of raw EMG
% 
% 
% Notes:    Ensure that events are places in the c3d file at the start & end
%           of each interval. e.g. 4 events == 3 intervals. If we want to set
%           the time vector to be the middle interval (between events 2&3), then
%           interval4Time = 2.
% 
%           e.g. multipleEMGprocess('myfile.c3d', 'emgset', 2)
% 
%           The output file will contain the time column (from interval4Time)
%           and the time normalized AVERAGE emg data over all intervals
% 
% --------------------------------------------------------------------

function [eVecGlob, EMGVecGlob] = multipleEMGprocess(c3dFile, emgSetName, emgProcessTasks, interval4Time)

usage = 'Usage: [eVecGlob, EMGVecGlob] = multipleEMGprocess(c3dFile, emgSetName, emgProcessTasks, interval4Time)';

if nargin ~= 4
    disp(usage)
    return
end

C3DkeyEMG = getEvents(c3dFile);
C3DkeyEMG.interval4Time = interval4Time;
loadLabels;
if ischar(emgProcessTasks)
    emgProcessTasks = glab.(emgProcessTasks);
end


% look for 'plot' and change value to C3Dkey to change the plotting option
l = length(emgProcessTasks);
ind = getIndex(emgProcessTasks, 'plot');
if ind == 0
    emgProcessTasks{l+1} = 'plot';
    emgProcessTasks{l+2} = 'C3Dkey';
else
    emgProcessTasks{ind+1} = 'C3Dkey';
end
    

% Save EMG
[eVecGlob, EMGVecGlob] = batchEMGprocess(C3DkeyEMG, emgSetName, emgProcessTasks, 'AverageEMG');


