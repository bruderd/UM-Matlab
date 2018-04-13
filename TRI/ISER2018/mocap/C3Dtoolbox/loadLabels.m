% This file contains parameters used to
%   1) Configure the scripts to the lab
%   2) Prepare labels for plotting
% ----------------------------------------
% LAB: LabName
% ----------------------------------------


% =================================
% SECTION 1: LAB DEPENDANCY
% =================================


% 1a) Force Plate Naming Convention
% ---------------------------------------------------------------------
% Force plate numbers MUST be set up within the lab.
% The names of the channels of these force plates are defined here.
%       EXAMPLE
%       glab.FP.string = '%s%d%s';
%       glab.FP.prefix  = {'FP','FP','FP','FP','FP','FP'};
%       glab.FP.suffix  = {'Fx','Fy','Fz','Mx','My','Mz'};
%       
%       Moment about FP #3 in the X direction analog channel = 'FP3Mx'
%       Please note that the glab.FP.string IS CASE SENSITIVE
% ---------------------------------------------------------------------
glab.FP.string = '%s%d%s';      % Prefix, PlateNum, Suffix

%%%% WORKSTATION
glab.FP.prefix  = {'FP','FP','FP','FP','FP','FP'};
glab.FP.suffix  = {'Fx','Fy','Fz','Mx','My','Mz'};

glab.FP.verticalForceIndex = 3;




% 1b) Coordinate Vector Setup
% ---------------------------------------------------------------------
% This section sets transformation vectors that are crucial to the
% rigid body transformations needed between the Force Plate, Vicon,
% and the Model coordinate systems.
% 
% Example: FP coord system -> Model coord system
%          glab.dirVec.FPMODEL(3,:)  = [2 -3 -1];   
% 
%          MODEL X = FP Y
%          MODEL Y = FP -Z
%          MODEL Z = FP -X
% ---------------------------------------------------------------------

% X direction (Vicon) Gait -> transformation vectors
glab.dirVec.FPMODEL(1,:)  = [-1 -3 -2];   % FP coord system    -> Model coord system
glab.dirVec.VICMODEL(1,:) = [1 3 -2];     % Vicon coord system -> Model coord system

% -X direction (Vicon) Gait -> transformation vectors
glab.dirVec.FPMODEL(2,:)  = [1 -3 2];   % FP coord system    -> Model coord system
glab.dirVec.VICMODEL(2,:) = [-1 3 2];   % Vicon coord system -> Model coord system

% Y direction (Vicon) Gait -> transformation vectors
glab.dirVec.FPMODEL(3,:)  = [2 -3 -1];   % FP coord system    -> Model coord system
glab.dirVec.VICMODEL(3,:) = [2 3 1];     % Vicon coord system -> Model coord system

% -Y direction (Vicon) Gait -> transformation vectors
glab.dirVec.FPMODEL(4,:)  = [-2 -3 1];   % FP coord system    -> Model coord system
glab.dirVec.VICMODEL(4,:) = [-2 3 -1];   % Vicon coord system -> Model coord system




% 1c) Vicon / OpenSim Offset Marker
% ---------------------------------------------------------------------
% This section sets marker to be used as the offset point to
% line up the model to the walking platform in OpenSim
% ---------------------------------------------------------------------
glab.offsetMarker = 'SACR';


% 1d) Joint / Muscle Model used in OpenSim
% ---------------------------------------------------------------------
% Set this corresponding to the joint / muscle names in the OpenSim model
% Currently the joint label is only used to get create the coordinates
% file in getKinetics.m
% ---------------------------------------------------------------------
glab.jointModel = 'jointDelp23';          % Joint convension
glab.muscleModel = 'muscleDelp92';         % Muscle convension






% 1e) Directory to store toolbox outputs
% ---------------------------------------------------------------------
% Output images, mat files, and text files from toolbox function calls
% can be stored in this directory for future reference. Must have
% a backslash at the end. i.e. '.\MyDir\'. This functionality is
% toggled by glab.storeInfo (1 = on, 0 = off)
% ---------------------------------------------------------------------
glab.storeInfo = 1;
glab.infoDirectory = '.\GaitExtract\';





% =================================
% SECTION 2: FORCE PLATE PLOTTING
% (Used by getKinetics.m)
% =================================


% 2a) Force Plate Labels
% ---------------------------------------------------------------------
glab.name{1} = 'GRF';
glab.name{2} = 'CoP';
glab.name{3} = 'GRMo';
glab.name{4} = 'GRMx';



% Force Plate Directions
% (x dir+-, y dir+-, z dir+-)
% ---------------------------------------------------------------------
glab.dir{1} = {'Fore',     'Aft'};              % X
glab.dir{2} = {'Vertical', ''};                 % Y
glab.dir{3} = {'Lateral',  'Medial'};           % Z




% GRF Labels
% (label, units, directionSign+, directionSign-)
% ---------------------------------------------------------------------
glab.S{1} = [{'GRF X (right)', 'N'}, glab.dir{1}(1), glab.dir{1}(2)];
glab.S{2} = [{'GRF Y (right)', 'N'}, glab.dir{2}(1), glab.dir{2}(2)];
glab.S{3} = [{'GRF Z (right)', 'N'}, glab.dir{3}(1), glab.dir{3}(2)];
glab.S{4} = [{'GRF X (left)',  'N'}, glab.dir{1}(1), glab.dir{1}(2)];
glab.S{5} = [{'GRF Y (left)',  'N'}, glab.dir{2}(1), glab.dir{2}(2)];
glab.S{6} = [{'GRF Z (left)',  'N'}, glab.dir{3}(1), glab.dir{3}(2)];




% CoP Labels
% (label, units, directionSign+, directionSign-)
% ---------------------------------------------------------------------
glab.X{1} = [{'CoP X (right)', 'm'}, glab.dir{1}(1), glab.dir{1}(2)];
glab.X{2} = [{'CoP Y (right)', 'm'}, glab.dir{2}(1), glab.dir{2}(2)];
glab.X{3} = [{'CoP Z (right)', 'm'}, glab.dir{3}(1), glab.dir{3}(2)];
glab.X{4} = [{'CoP X (left)',  'm'}, glab.dir{1}(1), glab.dir{1}(2)];
glab.X{5} = [{'CoP Y (left)',  'm'}, glab.dir{2}(1), glab.dir{2}(2)];
glab.X{6} = [{'CoP Z (left)',  'm'}, glab.dir{3}(1), glab.dir{3}(2)];




% GRM Labels (Origin)
% (label, units, directionSign+, directionSign-)
% ---------------------------------------------------------------------
glab.Mo{1} = [{'GRMo X (right)', 'Nm'}, glab.dir{1}(1), glab.dir{1}(2)];
glab.Mo{2} = [{'GRMo Y (right)', 'Nm'}, glab.dir{2}(1), glab.dir{2}(2)];
glab.Mo{3} = [{'GRMo Z (right)', 'Nm'}, glab.dir{3}(1), glab.dir{3}(2)];
glab.Mo{4} = [{'GRMo X (left)',  'Nm'}, glab.dir{1}(1), glab.dir{1}(2)];
glab.Mo{5} = [{'GRMo Y (left)',  'Nm'}, glab.dir{2}(1), glab.dir{2}(2)];
glab.Mo{6} = [{'GRMo Z (left)',  'Nm'}, glab.dir{3}(1), glab.dir{3}(2)];




% GRM Labels (CoP)
% (label, units, directionSign+, directionSign-)
% ---------------------------------------------------------------------
glab.Mx{1} = [{'GRMx X (right)', 'Nm'}, glab.dir{1}(1), glab.dir{1}(2)];
glab.Mx{2} = [{'GRMx Y (right)', 'Nm'}, glab.dir{2}(1), glab.dir{2}(2)];
glab.Mx{3} = [{'GRMx Z (right)', 'Nm'}, glab.dir{3}(1), glab.dir{3}(2)];
glab.Mx{4} = [{'GRMx X (left)',  'Nm'}, glab.dir{1}(1), glab.dir{1}(2)];
glab.Mx{5} = [{'GRMx Y (left)',  'Nm'}, glab.dir{2}(1), glab.dir{2}(2)];
glab.Mx{6} = [{'GRMx Z (left)',  'Nm'}, glab.dir{3}(1), glab.dir{3}(2)];



% =================================
% SECTION 3: OPENSIM MODEL LABELS
% =================================


% Joint Labels
% (label, unitsAngle, unitsTorque)
% Skeletal model used here: Anderson, F. C., 1999. A Dynamic Optimisation
% Solution for a Complete Cycle of Normal Gait, Ph.D Thesis, 
% The University of Texas at Austin, Austin.
% 17 joints, 6 spatial direction coords = 23 degrees of freedom
% ---------------------------------------------------------------------
glab.jointDelp23{1} = {'pelvis_tx',         'm',   'N'};
glab.jointDelp23{2} = {'pelvis_ty',         'm',   'N'};
glab.jointDelp23{3} = {'pelvis_tz',         'm',   'N'};
glab.jointDelp23{4} = {'pelvis_list',       'deg', 'Nm'};
glab.jointDelp23{5} = {'pelvis_rotation',   'deg', 'Nm'};
glab.jointDelp23{6} = {'pelvis_tilt',       'deg', 'Nm'};
glab.jointDelp23{7} = {'hip_flexion_r',     'deg', 'Nm'};
glab.jointDelp23{8} = {'hip_adduction_r',   'deg', 'Nm'};
glab.jointDelp23{9} = {'hip_rotation_r',    'deg', 'Nm'};
glab.jointDelp23{10} = {'knee_angle_r',     'deg', 'Nm'};
glab.jointDelp23{11} = {'ankle_angle_r',    'deg', 'Nm'};
glab.jointDelp23{12} = {'subtalar_angle_r', 'deg', 'Nm'};
glab.jointDelp23{13} = {'mtp_angle_r',      'deg', 'Nm'};
glab.jointDelp23{14} = {'hip_flexion_l',    'deg', 'Nm'};
glab.jointDelp23{15} = {'hip_adduction_l',  'deg', 'Nm'};
glab.jointDelp23{16} = {'hip_rotation_l',   'deg', 'Nm'};
glab.jointDelp23{17} = {'knee_angle_l',     'deg', 'Nm'};
glab.jointDelp23{18} = {'ankle_angle_l',    'deg', 'Nm'};
glab.jointDelp23{19} = {'subtalar_angle_l', 'deg', 'Nm'};
glab.jointDelp23{20} = {'mtp_angle_l',      'deg', 'Nm'};
glab.jointDelp23{21} = {'lumbar_extension', 'deg', 'Nm'};
glab.jointDelp23{22} = {'lumbar_bending',   'deg', 'Nm'};
glab.jointDelp23{23} = {'lumbar_rotation',  'deg', 'Nm'};






% =================================
% SECTION 4: EMG LABELS
% =================================
% (label, name)
% Used by batchEMGprocess.m to know which EMG analog channels to extract
% from the C3D file
% ---------------------------------------------------------------------


% Test EMG set
% --------------------
glab.testEMG{1}  = {'VASMED' ,'Vastus Medialis'};
glab.testEMG{2}  = {'GAS'    ,'Medial Gastrocnemius'};
glab.testEMG{3}  = {'SOL'    ,'Soleus'};
glab.testEMG{4}  = {'HAMLAT' ,'Lateral Hamstring'};




% EMG Processing Tasks
% {ProcessTask1, Value1, ProcessTask2, Value2, ...}
% notes: refer to processEMG.m for more information or type help processEMG
% 	     the operations are executed in order of appearance.
% ---------------------------------------------------------------------
glab.EMGprocessTKE1 = { 'tke', [], ...              % TKE filter
                        'rect', [], ...             % Rectification
                        'plot', 0 };                % Plot   
                    
                    
                    



% =================================
% SECTION 5: MARKERSET LABELS
% =================================
% Used by getMarkers.m to know which markers to extract
% from the C3D file
% ---------------------------------------------------------------------

glab.markersStatic = {'RSHO', 'LSHO', 'C7', ...
                      'RASI', 'LASI', 'SACR', ...
                      'RTHAP', 'RTHAD', 'RTHLD', 'RLEPI', 'RMEPI', ...
                      'RTIAP', 'RTIAD', 'RMMAL', 'RLMAL', ...
                      'RHEEL', 'RP1MT', 'RP5MT', 'RTOE', ...
                      'LTHAP', 'LTHAD', 'LTHLD', 'LLEPI', 'LMEPI', ...
                      'LTIAP', 'LTIAD', 'LMMAL', 'LLMAL', ...
                      'LHEEL', 'LP1MT', 'LP5MT', 'LTOE'};

                  

           
glab.markersDynamic = {'RSHO', 'LSHO', 'C7', ...
                      'RASI', 'LASI', 'SACR', ...
                      'RTHAP', 'RTHAD', 'RTHLD', 'RLEPI',  ...
                      'RTIAP', 'RTIAD', 'RLMAL', ...
                      'RHEEL', 'RP1MT', 'RP5MT', 'RTOE', ...
                      'LTHAP', 'LTHAD', 'LTHLD', 'LLEPI', ...
                      'LTIAP', 'LTIAD', 'LLMAL', ...
                      'LHEEL', 'LP1MT', 'LP5MT', 'LTOE'};             
                   
                   
 
