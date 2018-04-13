% Retrieve Key Events from a gait cycle in c3d format
% Tim Dorn
% May 2009
% 
% ----------------------------------------------------------------------
% Usage: C3Dkey = getEvents(c3dFile, direction*, FP_order*)
% ----------------------------------------------------------------------
% 
% Inputs:   c3dFile = the name of the c3d file
% 
%           direction* = the direction of travel in VICON COORDINATES.
%                           (required for marker extraction)
%                1 = X   ;  -1 = -X 
%                2 = Y   ;  -2 = -Y
% 
%           FP_order* = Vector containing the order of force plate numbers
%                      stepped on during the trial. This parameter can be 
%                      neglected if the trial is NOT a dynamic general
%                      trial (i.e (static trial or EMG/kinematics only trial) but
%                      otherwise it must be included.
%                      e.g. Walking on FP3, FP2, FP1 --> FP_order = [3 2 1]
%                      e.g. Running on FP1, FP3 (FP2 is skipped due to the
%                      large stride length --> FP_order = [1 3]
% 
%               Note: direction and FP_order are required here if you will
%                     be processing kinetics using this toolbox
% 
%                         
%                        
% 
% Outputs:  C3Dkey is a structure which contains key information about the trial
% 
%               C3Dkey.dirVec.FPMODEL: direction vector from FP coords -> Model coords
%               C3Dkey.dirVec.VICMODEL: direction vector from Vicon coords -> Model coords
%               C3Dkey.dirVec.FPVICON: direction vector from FP coords -> Vicon coords
% 
%               C3Dkey.dirVec.MODELFP: direction vector from Model coords -> FP coords
%               C3Dkey.dirVec.MODELVIC: direction vector from Model coords -> Vicon coords
%               C3Dkey.dirVec.VICONFP: direction vector from Vicon coords -> FP coords
% 
%               C3Dkey.name: subject name
%               C3Dkey.markerSet: marker set used for the trial
%               C3Dkey.c3dFile: c3d file name
%               C3Dkey.c3dFile: c3d file path
%               C3Dkey.direction: direction of forward facing value
%               C3Dkey.aFreq: analog frequency
%               C3Dkey.vFreq: video (VICON) frequency
%               C3Dkey.r: video/analog frequency ratio
%               C3Dkey.mass: subject mass (SIMPLE TRIAL ONLY)
% 
%               C3Dkey.numFrames.uncroppedV: number of total video frames in the trial
%               C3Dkey.numFrames.uncroppedA: number of total analog frames in the trial
%               C3Dkey.numFrames.croppedV = number of cropped frames (video)
%               C3Dkey.numFrames.croppedA = number of cropped frames (analog)
% 
%               C3Dkey.event.txt: label of events
%               C3Dkey.event.times: times of events
%               C3Dkey.event.percent: percentage of cycle that events occur
%               C3Dkey.event.Vframe: video frames of events
%               C3Dkey.event.Aframe: analog frames of events
%               C3Dkey.event.Vframe0: video frames of events (starting at frame 1)
%               C3Dkey.event.Aframe0: analog frames of events (starting at frame 1)
%               C3Dkey.event.times0: times of events (starting at time 0)
% 
%               C3Dkey.interval.txt: txt interval of events
%               C3Dkey.interval.time: time interval of events
%               C3Dkey.interval.Vframe: video frame interval of events
%               C3Dkey.interval.Aframe: analog frame interval of events
%               C3Dkey.interval.time0: time interval of events (starting at time 0)
%               C3Dkey.interval.Vframe0: video frame interval of events (starting at frame 1)
%               C3Dkey.interval.Aframe0: analog frame interval of events (starting at frame 1)
% 
%               C3Dkey.sequence.frames: force plate frame sequence
%               C3Dkey.sequence.plates: force plate number sequence
%               C3Dkey.sequence.eventIndex: force plate event boundaries
%               C3Dkey.sequence.txt: force plate event text sequence
% 
%               C3Dkey.offset: offset in mm to put the model on the platform in OpenSim
%               C3Dkey.averageSpeed: average trial speed (m/s)
%               C3Dkey.FP_order: Force plate order (in terms of stepping #)
%               C3Dkey.FP_order_inv: Force plate order inverse
%               C3Dkey.trialType: trial type (either will be SIMPLE or GENERAL)
%               C3Dkey.numPlatesTotal = total number of force plates in the trial
%               C3Dkey.numPlatesUsed = number of force plates used for extraction
%               C3Dkey.stanceFramesV = (1,:) - Right leg, (2,:) - Left leg (vicon frames)
% 
%               C3Dkey.allowed.markers:   1 if markers can be extracted from this c3dfile
%               C3Dkey.allowed.kinetics:  1 if markers can be extracted from this c3dfile
%               C3Dkey.allowed.EMG:       1 if EMG can be extracted from this c3dfile
% 
% 
%       Time Vectors (of labelled event):
%       ---------------------------------
%       C3Dkey.timeVec.c3dAnalogFrame:  actual analog frame number (from c3d)
%       C3Dkey.timeVec.analogFrame:     analog frame number (starting at 1)
%       C3Dkey.timeVec.c3dVideoFrame:   actual video frame number (from c3d)
%       C3Dkey.timeVec.videoFrame:      video frame number (starting at 1)
%       C3Dkey.timeVec.Asec:            analog time (sec) -> starting at 0 sec
%       C3Dkey.timeVec.Vsec:            video time (sec) -> starting at 0 sec
%       C3Dkey.timeVec.Apercent:        analog percentage of labeled event
%       C3Dkey.timeVec.Vpercent:        video percentage of labeled event 
%
% 
% Notes
% -----
% 
% Events must be examined and labeled in Vicon before using this script.
% 
% ----------------------------------------------------------------------

function C3Dkey = getEvents(c3dFile, direction, FP_order)

C3Dkey.allowed.markers = 1;
C3Dkey.allowed.kinetics = 1;
C3Dkey.allowed.EMG = 1;

if nargin == 1,
    direction = 0;
    FP_order = [];
    C3Dkey.allowed.markers = 0;
    C3Dkey.allowed.kinetics = 0;
    
elseif nargin == 2,
    FP_order = [];
    C3Dkey.allowed.kinetics = 0;
elseif nargin ~= 3,
    fprintf('Usage: getEvents(c3dFile, direction, FP_order*)\n');
    clear C3Dkey
    return
end

if abs(direction > 2)
    error('Invalid direction scalar: type help getEvents for allowable direction values');
end


% Open C3D File
% -------------
loadLabels;
[c3dpathstr, c3dname] = fileparts(c3dFile);
C3Dkey.c3dFile = c3dname;
C3Dkey.c3dPath = c3dpathstr;
itf = c3dserver();
openc3d(itf, 0, C3Dkey.c3dFile);
fprintf('\nEXTRACTING EVENTS FOR: %s\n--------------------------------------------\n', c3dFile);



% Extract Misc Trial Information
% ------------------------------

try
    aIndex = itf.GetParameterIndex('SUBJECTS', 'NAMES');
    C3Dkey.name = itf.GetParameterValue(aIndex, 0);
    if isempty(C3Dkey.name)
        C3Dkey.name = input('No Subject Name Found. Enter Name Here: ', 's');
    %     C3Dkey.name = 'SubjectName';
    end

    Units_index = itf.GetParameterIndex('POINT', 'UNITS');
    C3Dkey.units = itf.GetParameterValue(Units_index, 0);
    if strcmp(C3Dkey.units, 'mm')
        C3Dkey.divide_to_meters = 1000;
        C3Dkey.divide_to_mm = 1;
    else
        C3Dkey.divide_to_meters = 1;
        C3Dkey.divide_to_mm = 0.001;
    end
    
    aIndex = itf.GetParameterIndex('SUBJECTS', 'MARKER_SETS');
    C3Dkey.markerSet = itf.GetParameterValue(aIndex, 0);
    if isempty(C3Dkey.markerSet)
        C3Dkey.markerSet = 'MarkerSet';
    end
    fprintf('Subject Name: %s\nMarkerset: %s\n', C3Dkey.name, C3Dkey.markerSet);
catch
    fprintf('WARNING: Cameras have been disabled for this trial\n');
end
    

% Capture Frequencies
C3Dkey.r = double(itf.GetAnalogVideoRatio);
C3Dkey.vFreq = double(itf.GetVideoFrameRate);
C3Dkey.aFreq = C3Dkey.r * C3Dkey.vFreq;
C3Dkey.numFrames.uncroppedV = nframes(itf);
C3Dkey.numFrames.uncroppedA = C3Dkey.numFrames.uncroppedV * C3Dkey.r - C3Dkey.r + 1;


    

% Extract Coordinate & Force Plate Information
% --------------------------------------------
aIndex = itf.GetParameterIndex('FORCE_PLATFORM', 'USED');
C3Dkey.numPlatesTotal = double(itf.GetParameterValue(aIndex, 0));
C3Dkey.numPlatesUsed = length(FP_order);

C3Dkey.FP_order = FP_order;
C3Dkey.FP_order_inv = invFP(C3Dkey.FP_order);

C3Dkey.direction = direction;
ind = find([1 -1 2 -2] == C3Dkey.direction);

if ~isempty(ind),
    C3Dkey.dirVec.FPMODEL  = glab.dirVec.FPMODEL(ind,:);
    C3Dkey.dirVec.VICMODEL = glab.dirVec.VICMODEL(ind,:);

    % FP coord system -> Vicon coord system
    C3Dkey.dirVec.FPVICON = toVector(toMatrix(invDirVec(C3Dkey.dirVec.VICMODEL))* ...
        toMatrix(C3Dkey.dirVec.FPMODEL));

    % Inverse transformation vectors
    C3Dkey.dirVec.MODELFP = invDirVec(C3Dkey.dirVec.FPMODEL);
    C3Dkey.dirVec.MODELVIC = invDirVec(C3Dkey.dirVec.VICMODEL);
    C3Dkey.dirVec.VICONFP = invDirVec(C3Dkey.dirVec.FPVICON);
end



% Extract Event Information
% -------------------------
noEvents = 0;
try
    aIndex = itf.GetParameterIndex('EVENT', 'USED');
    C3Dkey.event.nEvent = double(itf.GetParameterValue(aIndex, 0));
catch
    noEvents = 1;
end

% Detected no events
if C3Dkey.event.nEvent == 0 || noEvents,
    fprintf('WARNING: No events detected in the c3d file.\n')
    insertArtEvents = input('Do you want to insert artificial events at the start/end frame? [y/n]: ', 's');
    if insertArtEvents == 'y'
        txtRaw = {'GEN', 'GEN'};
        timeRaw = [itf.GetVideoFrame(0), itf.GetVideoFrame(1)]/C3Dkey.vFreq;
        C3Dkey.event.nEvent = 2;
    else
        return
    end
    
% Detected some events
else
    fprintf('Detected %d events\n', C3Dkey.event.nEvent);
    bIndex = itf.GetParameterIndex('EVENT', 'CONTEXTS');
    cIndex = itf.GetParameterIndex('EVENT', 'LABELS');
    dIndex = itf.GetParameterIndex('EVENT', 'TIMES');

    for i = 1:C3Dkey.event.nEvent
       txtRawtmp = [itf.GetParameterValue(bIndex, i-1),...
               itf.GetParameterValue(cIndex, i-1)];
       timeRaw(i) = double(itf.GetParameterValue(dIndex, 2*i-1));

       if     strmatch(upper(txtRawtmp),'RIGHTFOOT OFF')    txtRawtmp= 'rFO'; 
       elseif strmatch(upper(txtRawtmp),'LEFTFOOT OFF')     txtRawtmp= 'lFO'; 
       elseif strmatch(upper(txtRawtmp),'RIGHTFOOT STRIKE') txtRawtmp= 'rFS';
       elseif strmatch(upper(txtRawtmp),'LEFTFOOT STRIKE')  txtRawtmp= 'lFS';
       elseif strmatch(upper(txtRawtmp),'RIGHTEVENT')       txtRawtmp= 'rGEN';
       elseif strmatch(upper(txtRawtmp),'LEFTEVENT')        txtRawtmp= 'lGEN';
       elseif strmatch(upper(txtRawtmp),'GENERALEVENT')     txtRawtmp= 'GEN';
       end

       txtRaw{i} = txtRawtmp;
    end
end


[timeNew, idNew] = sort(timeRaw);       % sort the events in time order
for i = 1:C3Dkey.event.nEvent;
    j = idNew(i);
    C3Dkey.event.txt{i} = txtRaw{j};
    time(i) = timeRaw(j);
end




% Create the list of event frame numbers
% --------------------------------------
% Note that a video frame V is equivalent to an analog frame A where
% A = (V-1)*r + 1, where r is the analog/video frame ratio. For example
% if r = 9 and we start at video frame 8, this is equivalent to analog
% frame (8-1)*9 + 1 = 64
% ----------------------------------------------------------------------

C3Dkey.event.times = time;
C3Dkey.event.times0 = C3Dkey.event.times - C3Dkey.event.times(1);

for i = 1:C3Dkey.event.nEvent,
   C3Dkey.event.percent(i) = (time(i) - time(1)) / (time(end) - time(1)) * 100;
end

C3Dkey.event.Vframes = round(C3Dkey.vFreq * time(1:C3Dkey.event.nEvent) + 1)';
C3Dkey.event.Aframes = C3Dkey.event.Vframes * C3Dkey.r - C3Dkey.r + 1;

C3Dkey.event.Vframes0 = C3Dkey.event.Vframes - C3Dkey.event.Vframes(1) + 1;
C3Dkey.event.Aframes0 = C3Dkey.event.Aframes - C3Dkey.event.Aframes(1) + 1;

C3Dkey.numFrames.croppedV = C3Dkey.event.Vframes(end) - C3Dkey.event.Vframes(1) + 1;
C3Dkey.numFrames.croppedA = C3Dkey.event.Aframes(end) - C3Dkey.event.Aframes(1) + 1;



% Get trial intervals
% -------------------
C3Dkey.interval.numInterval = C3Dkey.event.nEvent-1;
for i = 1:C3Dkey.event.nEvent-1;
    C3Dkey.interval.txt(i,:) = [C3Dkey.event.txt(i),C3Dkey.event.txt(i+1)];
    C3Dkey.interval.time(i,:) = [C3Dkey.event.times(i), C3Dkey.event.times(i+1)];
    C3Dkey.interval.Vframe(i,:) = [C3Dkey.event.Vframes(i), C3Dkey.event.Vframes(i+1)];
    C3Dkey.interval.Aframe(i,:) = [C3Dkey.event.Aframes(i), C3Dkey.event.Aframes(i+1)];
    C3Dkey.interval.time0(i,:) = [C3Dkey.event.times0(i), C3Dkey.event.times0(i+1)];
    C3Dkey.interval.Vframe0(i,:) = [C3Dkey.event.Vframes0(i), C3Dkey.event.Vframes0(i+1)];
    C3Dkey.interval.Aframe0(i,:) = [C3Dkey.event.Aframes0(i), C3Dkey.event.Aframes0(i+1)];
end



% Set up various time vectors
% ---------------------------
C3Dkey.timeVec.c3dAnalogFrame = C3Dkey.event.Aframes(1) : C3Dkey.event.Aframes(end);
C3Dkey.timeVec.analogFrame = C3Dkey.timeVec.c3dAnalogFrame - C3Dkey.event.Aframes(1) + 1;
C3Dkey.timeVec.c3dVideoFrame = C3Dkey.event.Vframes(1) : C3Dkey.event.Vframes(end);
C3Dkey.timeVec.videoFrame = C3Dkey.timeVec.c3dVideoFrame - C3Dkey.event.Vframes(1) + 1;
C3Dkey.timeVec.Asec = (C3Dkey.timeVec.analogFrame - 1) / C3Dkey.aFreq;
C3Dkey.timeVec.Vsec = (C3Dkey.timeVec.videoFrame - 1) / C3Dkey.vFreq;
C3Dkey.timeVec.Apercent = makepercent(C3Dkey.timeVec.analogFrame-1);
C3Dkey.timeVec.Vpercent = makepercent(C3Dkey.timeVec.videoFrame-1);



% Detect the offset to the platform in mm (OpenSim)
% and get the average trial speed (m/s)
% -------------------------------------------------
if strcmp(glab.offsetMarker, '');
    fprintf('No marker offset provided. Offset set to [0 0 0]\n')
    C3Dkey.offset = [0 0 0];
elseif C3Dkey.allowed.markers == 0
    fprintf('Warning: No direction given. Unable to compute marker offset.\n')
    fprintf('Marker offset set to [0 0 0]\n')
    C3Dkey.offset = [0 0 0];
else
    try
        tmp = -get3dtarget(itf, glab.offsetMarker, 0, C3Dkey.event.Vframes(1), C3Dkey.event.Vframes(end));
        C3Dkey.offset = coordChange(tmp(1,:)', C3Dkey.dirVec.VICMODEL)';
        C3Dkey.offset(2) = 0;
        fprintf('\n[%s Offset Calculation (%s)] --> [X=%.2f  Y=%.2f  Z=%.2f] mm\n', ...
            c3dFile, glab.offsetMarker, C3Dkey.offset(1), C3Dkey.offset(2), C3Dkey.offset(3));
        
        tmp2 = coordChange(tmp', C3Dkey.dirVec.VICMODEL)';
        C3Dkey.averageSpeed = double(abs((tmp2(end,1) - tmp2(1,1)) / 1000 / ...
            (C3Dkey.timeVec.Vsec(end) - (C3Dkey.timeVec.Vsec(1)))));
        fprintf('Average Trial Speed: %.3f m/s\n', C3Dkey.averageSpeed)
        tspeed = C3Dkey.averageSpeed;
        save trialSpeed.data tspeed -ASCII      % Write average speed to file
        
    catch
        C3Dkey.offset = [0 0 0];
        fprintf('\nERROR: Marker %s does not exist. Offset has been ignored...\n', ...
            glab.offsetMarker)      
    end   
end





% ==================================================
% Detect feet on force plates during event intervals
% Used by getKinetics.m for kinetic extraction
% ==================================================
C3Dkey.sequence.frames = C3Dkey.interval.Aframe - C3Dkey.event.Aframes(1) + 1;   %analog frames starting at 1
numSeq = length(C3Dkey.sequence.frames(:,1));
C3Dkey.sequence.plates = zeros(numSeq, 2);
C3Dkey.sequence.txt = C3Dkey.interval.txt;


if length(C3Dkey.event.txt) == 2
    
    % TWO LABELLED EVENTS
    % --------------------
    % Only start and end events are labelled.
    % This type of trial is used for static trials or trials where only
    % EMG has been recorded. In cases like this, we would like to try to
    % obtain the subject's mass from the force plate (if indeed it is a
    % static trial. 
    % -----------------------------------------------------------------
    C3Dkey.trialType = '2_EVENTS';      % Simple / Static Trial
    
    isStatic = input('2 Events are detected. Is this a static trial? [y/n]: ', 's');
    if isStatic == 'y'
        
        fprintf('Defined as a STATIC TRIAL\n');
        % Get static mass of subject (useful for static trials)
        [C3Dkey.mass, plateNum] = getStaticMass(itf, C3Dkey);
        C3Dkey.mass = double(C3Dkey.mass);
        C3Dkey.sequence.plates = [plateNum 0];
        fprintf('Subject Mass: %f kg\n', C3Dkey.mass);
        tmass = C3Dkey.mass;
        save tmass.data tmass -ASCII      % Write mass to file

        % Save zero coordinates file (for static trial in OpenSim)
        colnamesJoints = {'time'};
        datarows = length(C3Dkey.timeVec.Vsec');
        numJoints = length(glab.(glab.jointModel));

        file = sprintf('%s_coordinates.mot', C3Dkey.c3dFile);
        bigM = [C3Dkey.timeVec.Vsec', zeros(datarows, numJoints)];
        for i = 1:numJoints
            colnamesJoints(i+1) = glab.(glab.jointModel){i}(1);
        end
        generateMotFile(bigM, colnamesJoints, file);
    else
        fprintf('Defined as a NON-STATIC TRIAL\n');
        C3Dkey = detectSeq(C3Dkey, itf, numSeq);
    end

else
      
    % GENERAL TRIAL PATTERN WITH > 2 LABELLED EVENTS
    % ----------------------------------------------
    % This could be a trial where start, end and intermediate events are labelled.
    % Alternativly, it could identify the boundaries of multiple gait cycles.
    % Either way, we want to try and detect which foot lands on which plate.
    % -----------------------------------------------------------------
    C3Dkey.trialType = '>2_EVENTS';           % General / Dynamic Trial
    
    if C3Dkey.numPlatesUsed > C3Dkey.numPlatesTotal,
        error('Error in FP_order: number of plates in FP_order exceeds number of plates in the trial...')
    end
    
    if isempty(FP_order)
       fprintf('Warning: No force plate order vector given (FP_order). Skipping FP detection.\n')
       fprintf('Warning: Can''t use this event key to extract kinetics...\n')
        
    else
    
        % Algorithm to detect left and right foot plate force plate orders.
        % To determine if a foot hits a plate, we check the vertical FP 
        % averages across the plates for each sequence, and check that they
        % are above a cricical minimum value (vertForceCutoff).
        C3Dkey = detectSeq(C3Dkey, itf, numSeq);
    end
end 

fprintf('Trial %s detected as: %s\n', c3dFile, upper(C3Dkey.trialType));

tmpallowed = [C3Dkey.allowed.markers, C3Dkey.allowed.kinetics, C3Dkey.allowed.EMG];
tmpallowedTxt1 = {'MARKERS', 'KINETICS', 'EMG'};
tmpallowedTxt2 = {'NO', 'YES'};
fprintf('==============================================\n');
for i = 1:length(tmpallowed)
    fprintf('Key compatible for extracting %s: %s\n', ...
        tmpallowedTxt1{i}, tmpallowedTxt2{tmpallowed(i)+1});
end
fprintf('==============================================\n');

closec3d(itf);


% Save event information to file
% ------------------------------
if glab.storeInfo == 1,
    if exist(glab.infoDirectory, 'dir') ~= 7,
        mkdir(glab.infoDirectory);
    end
    
    f1 = sprintf('%sKey_%s.mat', glab.infoDirectory, C3Dkey.c3dFile);
    save(f1, 'C3Dkey')
    
    % save txt file of events...
    filename = sprintf('%sEVENTS_%s.txt', glab.infoDirectory, C3Dkey.c3dFile);
    fid = fopen(filename, 'w');
    if fid < 0
        fprintf('\nERROR: %s could not be opened for writing...\n\n', filename);
        return
    end
    
    fprintf(fid, 'EVENT        time           time0         Percent         Vframes        Vframes0         Aframes        Aframes0\n');
    for i = 1:length(C3Dkey.event.txt)
       fprintf(fid, '%s\t\t%9.3f\t\t%9.3f\t\t%9.3f\t\t%9.3f\t\t%9.3f\t\t%9.3f\t\t%9.3f\n', ...
           C3Dkey.event.txt{i}, C3Dkey.event.times(i), C3Dkey.event.times0(i), ...
           C3Dkey.event.percent(i), C3Dkey.event.Vframes(i), C3Dkey.event.Vframes0(i), ...
           C3Dkey.event.Aframes(i), C3Dkey.event.Aframes0(i));     
    end
    fclose(fid);
end




% ========================================================================
% SUBFUNCTION: C3Dkey = detectSeq(C3Dkey, itf, numSeq)
% ========================================================================
% Detect foot sequence in C3Dkey. Adds numbers to the C3Dkey.stanceFramesV
% and C3Dkey.sequence fields
% 
% for each sequence in the trial,
% work out which plates have force on them (these are the active plates)
% assign the active plates to the forcePlateOrder
% if no plates are active, then we are airborne
% ========================================================================
    

function C3Dkey = detectSeq(C3Dkey, itf, numSeq)

vertForceCutoff = 10;    % Cutoff to detect if foot is on/off a plate (N)
R=1; L=2;


try
for i = 1:numSeq

    % determine active plate numbers in the sequence
    for j = 1:C3Dkey.numPlatesUsed
        avrF(j) = avrVertForceOnPlate(itf, abs(C3Dkey.FP_order(j)), ...
            C3Dkey.interval.Vframe(i,1), C3Dkey.interval.Vframe(i,2));
    end
    activePlates = abs(C3Dkey.FP_order(find(avrF>vertForceCutoff)));

    % zero support sequence (airborne)
    if isempty(activePlates),
        C3Dkey.sequence.plates(i,:) = [0 0];

    % single support sequence
    elseif length(activePlates) == 1,
        if findstr(C3Dkey.sequence.txt{i,1}, 'FO')
            footOnPlate = opp(C3Dkey.sequence.txt{i}(1));
        elseif findstr(C3Dkey.sequence.txt{i,1}, 'FS')
            footOnPlate = C3Dkey.sequence.txt{i,1}(1);
        elseif findstr(C3Dkey.sequence.txt{i,1}, 'GEN')
            footOnPlate = C3Dkey.sequence.txt{i,1}(1);
        else
            error('unable to detect single foot activePlate...')
        end

        switch footOnPlate
            case 'r'
                C3Dkey.sequence.plates(i,R) = activePlates;
            case 'l'
                C3Dkey.sequence.plates(i,L) = activePlates;
        end

    % double support sequence    
    else

        % find foot-off foot
        for k = 1:2
            if findstr(C3Dkey.sequence.txt{i,k}, 'FO')
                FOfoot = C3Dkey.sequence.txt{i,k}(1);
            end
        end

        % sort activePlates into gait style order
        if find(abs(C3Dkey.FP_order)==activePlates(1)) > ...
                find(abs(C3Dkey.FP_order)==activePlates(2))
            activePlatesOrdered = fliplr(activePlates);
        else
            activePlatesOrdered = activePlates;
        end

        % the foot representing FS is on the front force plate and
        % the foot representing FO is on the back force plate 
        if strcmp(FOfoot, 'r')
            % right foot at back
            C3Dkey.sequence.plates(i,:) = activePlatesOrdered;
        else
            % left foot at back
            C3Dkey.sequence.plates(i,:) = fliplr(activePlatesOrdered);
        end
    end
end
%     disp(C3Dkey.sequence.plates)     % DEBUG STEP


% Use C3Dkey.sequence.plates to determine the event boundaries
% for each force plate used in the trial
% -------------------------------------------------------------
for i = 1:C3Dkey.numPlatesUsed
    [ii,jj] = find(C3Dkey.sequence.plates == C3Dkey.FP_order(i));
    C3Dkey.sequence.eventIndex(i,:) = [ii(1) ii(end)+1];
end


catch
    fprintf('Warning: Could not detect force plate sequence.\n');
end






% Detect Stance Frames
% ---------------------
try
    C3Dkey.stanceFramesV.right = ...
        getTimes('rFS', 'rFO', C3Dkey.event.txt, C3Dkey.event.Vframes);
catch
    fprintf('Warning: Could not detect right foot stance frames.\n');
end

try
    C3Dkey.stanceFramesV.left = ...
        getTimes('lFS', 'lFO', C3Dkey.event.txt, C3Dkey.event.Vframes);  
catch
    fprintf('Warning: Could not detect left foot stance frames.\n');
end












% ========================================================================
% SUBFUNCTION: avrF = avrVertForceOnPlate(itf, plateNum, startFrame, endFrame)
% ========================================================================
% Returns the average vertical force of a force plate over a specified
% time period.

function avrF = avrVertForceOnPlate(itf, plateNum, startVFrame, endVFrame)

loadLabels;
plateLabel = sprintf(glab.FP.string, glab.FP.prefix{glab.FP.verticalForceIndex}, ...
    plateNum, glab.FP.suffix{glab.FP.verticalForceIndex});
try
    warning off MATLAB:divideByZero
    forceVec = getanalogchannel(itf, plateLabel, startVFrame, endVFrame);
%     if ~isempty(forceVec),
        avrF = mean(abs(forceVec));
%     end
catch
    fprintf('ERROR: Force Plate [%s] does not exist.\n', plateLabel);
    return     
end

% DEBUG STEP
% figure;
% plot(abs(forceVec));
% title(sprintf('plate %d frame %d:%d', plateNum, startVFrame, endVFrame));
% hline(avrF);





% ========================================================================
% SUBFUNCTION: m = toMatrix(dirVec)
% ========================================================================
% Turns the rigid body direction vector into a rigid body transformation 
% matrix

function m = toMatrix(dirVec)

m = zeros(3,3);
for i = 1:3
    m(i,abs(dirVec(i))) = sign(dirVec(i)); 
end





% ========================================================================
% SUBFUNCTION: m = toVector(m)
% ========================================================================
% Turns the rigid body transformation matrix into a rigid body direction
% vector 

function dirVec = toVector(m)

dirVec = zeros(1,3);
for i = 1:3
    for j = 1:3
        if abs(m(i,j)) == 1
            pos = j;
            sgn = sign(m(i,j));
        end
    end
    dirVec(i) = pos*sgn; 
end





% ========================================================================
% SUBFUNCTION: invDirVec
% ========================================================================
% Returns the inverse of a rigid body transformation vector

function dirVecInv = invDirVec(dirVec)

dirVecInv = toVector(inv(toMatrix(dirVec)));





% ========================================================================
% SUBFUNCTION: invFP
% ========================================================================
% Returns the inverse of the given FP_order vector
% Input: FP_order (i.e position 1 = first vicon plate number stepped on)
% Output: FP_inv  (i.e. position 1 = order of stepping of vicon plate 1)

function fpinv = invFP(fp)

if isempty(fp),
    fpinv = [];
else
    for i = 1:length(fp)
        fpinv(abs(fp(i))) = i;
    end
end





% ========================================================================
% SUBFUNCTION: mass = getStaticMass(itf, C3Dkey)
% ========================================================================
% Gets the mass of the subject. This should be only used during static
% trials for reasonable results. This function examines all channels of all
% force plates between the event periods and takes an average of the force
% on each channel. The maximum average would correspond to the subject's
% downward vertical force and hence divinging by 9.81 represents its mass

function [mass, plateNum] = getStaticMass(itf, C3Dkey)

g = 9.81;
for p = 1:C3Dkey.numPlatesTotal
  try
    avrF(p) = avrVertForceOnPlate(itf, p, ...
        C3Dkey.event.Vframes(1), C3Dkey.event.Vframes(end));
  catch
    fprintf('ERROR: Force plate #%d does not exist on this trial\n', p)
     mass = input('Enter subject mass in kg (press ENTER to ignore mass): ');
     if ~isnumeric(mass),
         mass = 0;
         fprintf('Invalid Mass: Mass will be set to 0. You can edit this in C3Dkey.mass later\n')
     end
     plateNum = -1;
     return
  end
end
[force, plateNum] = max(avrF);
mass = force / g;

detectionMass = 10;     % kg
if mass < detectionMass,
    fprintf('Warning: Detected subject mass (%.3f kg) is < %d kg\n', mass, detectionMass)
    newmass = input('Enter new subject mass in kg or press ENTER to accept the current mass: ');
    if ~isnumeric(newmass),
         mass = 0;
         plateNum = -1;
         fprintf('Invalid Mass: Mass will be set to 0. You can edit this in C3Dkey.mass later\n')
    elseif isempty(newmass),
        return
    else mass = newmass;
    end
end   
 




% ========================================================================
% SUBFUNCTION: times = getTimes(from, to, events, frames)
% ========================================================================
% Gets the frame values from and to a particular event.

function times = getTimes(from, to, events, frames)

if length(frames) ~= length(events)
    error('Event and frame vector must be the same length...')
end

% find the index of the first occuring 'from' event
[fromNum, fromIndex] = getFrame(from, events, frames);
fromIndex = fromIndex(1);
fromNum = fromNum(1);

% find the index of the first 'to' event that follows after the 'from' event
[toNum, toIndex] = getFrame(to, events, frames);
for i = 1:length(toNum)
    if toIndex(i) > fromIndex,
        toNum = toNum(i);
        toIndex = toIndex(i);
        times = [fromNum, toNum];
        return
    end
end

% if we get here, the 'to' event actually comes before the 'from' event
% in which case we need to wrap and have 2 parts to the frame values
toIndex = toIndex(1);
toNum = toNum(1);
times = [fromNum, frames(end) ; frames(1), toNum];





% ========================================================================
% SUBFUNCTION: [frameNum, frameIndex] = getFrame(event, eventList, frameList)
% ========================================================================
% Gets the frame number & index of a particular event from a frameList

function [frameNum, frameIndex] = getFrame(event, eventList, frameList)

l = length(eventList);
frameIndex = [];
frameNum = [];

for i = 1:l,
    if strcmp(event, eventList{i})
        frameIndex = [frameIndex, i];
        frameNum = frameList(frameIndex);
    end
end





% ========================================================================
% SUBFUNCTION: opp
% ========================================================================
% Used to detect the gait style.

function p2 = opp(p1)

if p1 == 'l', p2 = 'r'; end
if p1 == 'r', p2 = 'l'; end




% ========================================================================
% SUBFUNCTION: strideLength = getStride(itf, C3Dkey, RfootMarker, LfootMarker)
% ========================================================================
% Gets the stride length for each leg using the HEEL markers
% Not used anymore...

% function strideLength = getStride(itf, C3Dkey)
% 
% R = 1;
% L = 2;
% 
% fprintf('Calculating Stride Length using HEEL markers\n');
% heelMarkers = [coordChange(get3dtarget(itf, RfootMarker, 0, C3Dkey.event.Vframes(1), C3Dkey.event.Vframes(end))', C3Dkey.dirVec.VICMODEL)', ...
%                coordChange(get3dtarget(itf, LfootMarker, 0, C3Dkey.event.Vframes(1), C3Dkey.event.Vframes(end))', C3Dkey.dirVec.VICMODEL)'];
%            
% heelMarkersForeAft = [heelMarkers(:,1)/C3Dkey.divide_to_meters, heelMarkers(:,4)/C3Dkey.divide_to_meters];
% 
% rFS = getTimes('rFS', 'lFS', C3Dkey.event.txt, C3Dkey.event.Vframes0);
% lFS = getTimes('lFS', 'rFS', C3Dkey.event.txt, C3Dkey.event.Vframes0);
% strideLength(R) = abs(heelMarkersForeAft(rFS(1), R) - heelMarkersForeAft(rFS(2), L));
% strideLength(L) = abs(heelMarkersForeAft(lFS(1), L) -
% heelMarkersForeAft(lFS(2), R));

