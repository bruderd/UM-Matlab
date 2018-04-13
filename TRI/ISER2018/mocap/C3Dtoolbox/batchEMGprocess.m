% Batch Process EMG signals
% Tim Dorn
% November 2008
% 
% --------------------------------------------------------------------
% Usage: [eVecGlob, EMGVecGlob] = batchEMGprocess(C3Dkey, emgSetName, emgProcessTasks, fileSuffix*)
% --------------------------------------------------------------------
% 
% Inputs:   C3Dkey = key of dynamic C3D file (from getEvents.m)
% 
%           emgSetName: the label of EMG names contained in the EMG set
%               (this must be defined in loadlabels.m as a cell: glab.[emgSetName])
% 
%           emgProcessTasks: the EMG processing options in the order of execution
%               (this must be defined in loadlabels.m as a cell: glab.[emgProcessTasks])
% 
%           fileSuffix* = suffix of mot file that is saved    
%                   if this is not included, or empty, then file is not saved
% 
% 
% Outputs:  eVecGlob = structure of processed EMG
%           EMGVecGlob = structure of raw EMG
% 
% 
% Notes:    N/A
% 
% --------------------------------------------------------------------

function [eVecGlob, EMGVecGlob] = batchEMGprocess(C3Dkey, emgSetName, emgProcessTasks, fileSuffix)

usage = 'Usage: [eVecGlob, EMGVecGlob] = batchEMGprocess(C3Dkey, emgSetName, emgProcessTasks, fileSuffix*)';
    
if nargin == 3,
    fileSuffix = [];

elseif nargin ~= 4,
    disp(usage)
    return
end

loadLabels;
EMGset = glab.(emgSetName);
for i = 1:length(EMGset)
    musc_name_orig = EMGset{i}{1};
    musc_name = removeSpaces(strtok(musc_name_orig, '.()'));
end

if ischar(emgProcessTasks)
    emgProcessTasks = glab.(emgProcessTasks);
end

ind = getIndex(emgProcessTasks, 'plot');
if strcmpi(emgProcessTasks{ind+1}, 'C3Dkey')
    emgProcessTasks{ind+1} = C3Dkey;
end
ind = getIndex(emgProcessTasks, 'vertlines');
if strcmpi(emgProcessTasks{ind+1}, 'C3Dkey')
    emgProcessTasks{ind+1} = C3Dkey;
end



% ---------------
% Fixed Variables
% ---------------

% Get Dynamic EMG Values
bigM = [];
p = 0;
colnames = {'time'};
EMGVecGlob = extractRawEMG(C3Dkey, emgSetName);
names = fieldnames(EMGVecGlob);
numMuscles = length(names);
file = sprintf('%s_%s.mot', C3Dkey.c3dFile, fileSuffix);



% Process all EMG channels
% ------------------------
for i = 1:numMuscles
    l = length(glab.(emgSetName){i});
    musc_name = names{i};
    
    if l < 3
        nM = 1;
        musc_name_in_model = {musc_name};
    else
        nM = l-2;
        musc_name_in_model = glab.(emgSetName){i}(3:3+nM-1);
    end
    
    EMGVecTmp = EMGVecGlob.(musc_name);
    
    for j = 1:nM
        
        % Process EMGVec(t) -> Processed EMG e(t)
        % ---------------------------------------
        eVecTmp = processEMG(EMGVecTmp, emgProcessTasks);


        % Set variables
        eVecGlob.(musc_name) = eVecTmp;
        if isfield(eVecTmp, 'SplitStrides')     % multiple cycles in a single C3Dfile
                                                % uses AVERAGE values here
            colnames = [colnames, {sprintf('%s.Raw', musc_name_in_model{j})}, ...
                {sprintf('%s.Processed', musc_name_in_model{j})}];
            bigM = [bigM, eVecGlob.(musc_name).SplitStrides.rawAvr, eVecGlob.(musc_name).SplitStrides.processedAvr];
            timeAvr = eVecGlob.(musc_name).SplitStrides.timeAvr1;
            p = 1;
            
        else        % normal use: single cycle in a single C3Dfile
            colnames = [colnames, {sprintf('%s.Raw', musc_name_in_model{j})}, ...
                {sprintf('%s.Processed', musc_name_in_model{j})}];
            bigM = [bigM, EMGVecGlob.(musc_name).data, eVecGlob.(musc_name).data];
        end
    end
end


if ~isempty(fileSuffix)
    if p == 0
        bigM = [C3Dkey.timeVec.Asec', bigM];       % add time column
        fprintf('Saving EMG values\n')
    else
        bigM = [timeAvr, bigM];       % add time column
        fprintf('Saving AVERAGE EMG values\n')
    end
    generateMotFile(bigM, colnames, file);
end

