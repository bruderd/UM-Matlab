% Write XML files for use in OpenSim.
% Tim Dorn
% Nov 2008
% 
% --------------------------------------------------------------------
% Usage: xmlString = writeXML(type, C3Dkey*)
% --------------------------------------------------------------------
% 
% Inputs:  type: the type of XML setup file to be created (case sensitive)
%               type = 'scale' --> create scale XML file
%               type = 'ik' --> create inverse kinematics XML setup file
%               type = 'id' --> create inverse dynamics XML setup file
%               type = 'rra' --> create residual reduction (pass 1) XML setup file
%               type = 'cmc' --> create computed muscle control XML setup file
%               type = 'forward' --> create forward dynamics XML setup file
%
%           C3Dkey*: the C3D key structure from getEvents
%                   (if not given, default parameters are used. These
%                    can be then modified manually using an XML editor)
% 
%           saveFile: 0 = don't save file, 1 = save file (default)
% 
%           
% 
% Outputs:  xmlString: the generated XML string
% 
%           the XML file is saved in the current working directory as
%           [C3Dkey.c3dFile]_Setup_[type].xml
% 
% 
% Notes:    The output files are templates only. They may need to be
%           altered to conform to the precise simulation (i.e. paths, 
%           filtering frequencies, input files, etc.
% 
% --------------------------------------------------------------------


function xmlString = writeXML(type, C3Dkey, saveFile)

usage = 'Usage: xmlString = writeXML(type, C3Dkey*)';
types = 'Valid Types: ''scale'', ''ik'', ''id'', ''rra'', ''cmc'', ''forward''';

if nargin == 1,
    C3Dkey = [];
    saveFile = 1;
    
elseif nargin == 2,
    saveFile = 1;
   
elseif nargin ~= 3,
    fprintf('%s\n%s\n', usage, types);
    return
end


if isempty(C3Dkey)
    % use default C3Dkey for manual editing
    C3Dkey.c3dFile = 'C3DFILE';
    C3Dkey.name = 'SUBJECTNAME';
    C3Dkey.timeVec.Asec = [0 9999];
end



% =============================
% Simulation XML Setting Files
% =============================
modelPrefix = '..\MU2392_arms';


% Model Files
% -----------
F.unscaledModelFile = sprintf('%s_1xStrength_v18.osim', modelPrefix);
F.scaledModelFile = sprintf('..\\%s_SCALED.osim', C3Dkey.name);
F.scaledAdjustedModelFile = sprintf('%s_SCALED_ADJUSTED.osim', C3Dkey.name);


% Experimental / Model Output Data Files
% --------------------------------------

F.C3DFILE = C3Dkey.c3dFile;      % c3d filename without extension
F.manualScaleSetFile = sprintf('%s_Scale_ScaleSet.xml', F.C3DFILE);
F.markerFile = sprintf('%s.trc', F.C3DFILE);
F.coordinateFile = sprintf('%s_coordinates.mot', F.C3DFILE);
F.kineticsFile = sprintf('%s_kinetics.mot', F.C3DFILE);
F.IKFile = sprintf('%s_ik_arms.mot', F.C3DFILE);
F.controlsFile = sprintf('%s_controls.xml', F.C3DFILE);
F.statesFile = sprintf('%s_states.sto', F.C3DFILE);



% Model & Simulation XML Files
% ----------------------------
F.markerSetFile = sprintf('%s_Scale_MarkerSet.xml', modelPrefix);
F.measurementSetFile = sprintf('%s_Scale_MeasurementSet.xml', modelPrefix);

F.taskSetFile.SCALE = sprintf('%s_Scale_Tasks.xml', modelPrefix);
F.taskSetFile.IK = sprintf('%s_IK_Tasks.xml', modelPrefix);
F.taskSetFile.RRA = sprintf('%s_RRA_Tasks.xml', modelPrefix);
F.taskSetFile.CMC = sprintf('%s_CMC_Tasks.xml', modelPrefix);

F.actuatorSetFile.RRA = sprintf('%s_RRA_Actuators.xml', modelPrefix);
F.actuatorSetFile.CMC = sprintf('%s_CMC_Actuators.xml', modelPrefix);

F.controlConstraintsFile.RRA = sprintf('%s_RRA_ControlConstraints.xml', modelPrefix);
F.controlConstraintsFile.CMC = sprintf('%s_CMC_ControlConstraints.xml', modelPrefix);


% Result Directories
% ------------------
F.ID_Results = '.\INVDYN\INVDYN_Results';
F.RRA_Results = 'RRA_Results';
F.CMC_Recults = 'CMC_Results';


% ================
% Common Settings
% ================

% Integrator Settings
% -------------------
S.inte.maximum_number_of_integrator_steps = 30000;
S.inte.maximum_integrator_step_size = 0.0001;
S.inte.integrator_error_tolerance = 5e-6;
S.inte.integrator_fine_tolerance = 5e-8;
S.inte.output_precision = 20;


% Optimizer Settings
% ------------------
S.optIK.optimizer_algorithm = 'cfsqp';      % ipopt or cfsqp

S.optRRA.use_fast_optimization_target = 'false';
S.optRRA.optimizer_derivative_dx = 0.0001;
S.optRRA.optimizer_convergence_criterion = 1e-6;
S.optRRA.optimizer_max_iterations = 2000;
S.optRRA.optimizer_print_level = 0;
S.optRRA.optimizer_algorithm = 'ipopt';      % ipopt or cfsqp

S.optCMC.use_fast_optimization_target = 'true';
S.optCMC.optimizer_derivative_dx = 1;
S.optCMC.optimizer_convergence_criterion = 1e-4;
S.optCMC.optimizer_max_iterations = 2000;
S.optCMC.optimizer_print_level = 0;
S.optCMC.optimizer_algorithm = 'ipopt';      % ipopt or cfsqp


% Other Settings
% --------------
S.kinematics_filter_freq = 8;        % Hz




% ====================================================================
% Based on the type of XML file to create, the data structure is built
% and saved to an xml file
% ====================================================================
timeRange = [C3Dkey.timeVec.Vsec(1), C3Dkey.timeVec.Vsec(end)];

switch type
    
    % SCALE XML FILE (SCALE)
    % ----------------------
    case 'scale',
        root = 'ScaleTool';
        
        % ScaleTool
        data.ATTRIBUTE.name = sprintf('%s_SCALED', C3Dkey.name);
        if isfield(C3Dkey, 'mass')
            data.mass = C3Dkey.mass;
        else
            data.mass = '9999999';
        end
        data.height = '9999999';
        data.age = '9999999';
        data.notes = 'enter notes';
        
        % ScaleTool -> GenericModelMaker
        data.GenericModelMaker.ATTRIBUTE.name = '';
        data.GenericModelMaker.model_file = F.unscaledModelFile;
        data.GenericModelMaker.marker_set_file = F.markerSetFile;
        
        % ScaleTool -> ModelScaler
        data.ModelScaler.ATTRIBUTE.name = '';
        data.ModelScaler.apply = 'true';
        data.ModelScaler.scaling_order = 'measurements manualScale';
        
%         the following line is used for manualScaling (if necessary)
%         data.ModelScaler.ScaleSet.ATTRIBUTE.file = F.manualScaleSetFile;

        data.ModelScaler.MeasurementSet.ATTRIBUTE.file = F.measurementSetFile;
        data.ModelScaler.marker_file = F.markerFile;
        data.ModelScaler.time_range = num2str(timeRange);
        data.ModelScaler.preserve_mass_distribution = 'true';
%         data.ModelScaler.output_model_file = F.scaledModelFile;         % output model without markers
        data.ModelScaler.output_scale_file = sprintf('%s_scaleSet_applied.xml', C3Dkey.name);
        
        % ScaleTool -> MarkerPlacer
        data.MarkerPlacer.ATTRIBUTE.name = '';
        data.MarkerPlacer.apply = 'true';
        data.MarkerPlacer.optimizer_algorithm = S.optIK.optimizer_algorithm;
        data.MarkerPlacer.IKTaskSet.ATTRIBUTE.file = F.taskSetFile.SCALE;
        data.MarkerPlacer.marker_file = F.markerFile;
        data.MarkerPlacer.coordinate_file = F.coordinateFile;
        data.MarkerPlacer.time_range = num2str(timeRange);
        data.MarkerPlacer.output_model_file = F.scaledModelFile;          % output model with markers
        data.MarkerPlacer.output_motion_file = sprintf('%s_static_output.mot', C3Dkey.name);
        
        
        
        
    % INVERSE KINEMATICS XML FILE (IK)
    % --------------------------------
    case 'ik',
        
        % IKTool
        root = 'IKTool';
        data.ATTRIBUTE.name = F.C3DFILE;
        data.model_file = F.scaledModelFile;
        data.optimizer_algorithm = S.optIK.optimizer_algorithm;
        data.IKTaskSet.ATTRIBUTE.file = F.taskSetFile.IK;
        
        % IKTool -> IKTrialSet
        data.IKTrialSet.ATTRIBUTE.name = '';
        data.IKTrialSet.objects.IKTrial.ATTRIBUTE.name = F.C3DFILE;
        data.IKTrialSet.objects.IKTrial.marker_file = F.markerFile;
        data.IKTrialSet.objects.IKTrial.coordinate_file = F.coordinateFile;
        data.IKTrialSet.objects.IKTrial.time_range = num2str(timeRange);
        data.IKTrialSet.objects.IKTrial.output_motion_file = F.IKFile;
        data.IKTrialSet.objects.IKTrial.include_markers = 'true';
    
     
        
        
    % INVERSE DYNAMICS XML FILE (ID)
    % ------------------------------
    case 'id',
        
        % AnalyzeTool
        root = 'AnalyzeTool';
        data.ATTRIBUTE.name = sprintf('%s_INVDYN_IK_arms', F.C3DFILE);
        data.model_file = F.scaledModelFile;
        data.replace_actuator_set = 'false';
        data.actuator_set_files = ' ';
        data.results_directory = F.ID_Results;
        data.output_precision = S.inte.output_precision;
        timeRange(1) = timeRange(1) + 0.001;
        data = addInitialFinalTimes(data, timeRange);
        data.coordinates_file = F.IKFile;
        data.lowpass_cutoff_frequency_for_coordinates = S.kinematics_filter_freq;
        data = addExternalLoadSettings(data, S, F);

        % AnalyzeTool -> AnalysisSet
        data.AnalysisSet.ATTRIBUTE.name = 'Analyses';
        data = addInverseDynamicsAnalysis(data);

  
        
        
    % RESIDUAL REDUCTION
    % ------------------       
    case 'rra',
    
        % CMCTool
        root = 'CMCTool';
        data.ATTRIBUTE.name = F.C3DFILE;
        
        % Settings Required for RRA Pass 1
        % (Determinate IVD with ideal actuators to calculate residuals 
        % & adjust the model accordingly)
        
        data.model_file = F.scaledModelFile;
        data.replace_actuator_set = 'true';
        data.actuator_set_files = F.actuatorSetFile.RRA;
        data.results_directory = F.RRA_Results;
        data = addInitialFinalTimes(data, timeRange);
%         data.initial_time_for_com_adjustment = timeRange(1);
%         data.final_time_for_com_adjustment = timeRange(end);
                
        data.compute_average_residuals = 'true';
        data.adjust_com_to_reduce_residuals = 'true';
        data.adjusted_com_body = 'torso';
        data.output_model_file = F.scaledAdjustedModelFile;        
        data.adjust_kinematics_to_reduce_residuals = 'true';
        
        data.desired_kinematics_file = F.IKFile;
        data.lowpass_cutoff_frequency = S.kinematics_filter_freq;
        data = addExternalLoadSettings(data, S, F);

        data = addIntegratorSettings(data, S.inte);

        
        
        % Settings Required for RRA Pass 2
        % (IVD Optimization to adjust actuators & kinematics)
        
        data.task_set_file = F.taskSetFile.RRA;
        data.constraints_file = F.controlConstraintsFile.RRA;       
        data.cmc_time_window = 0.001;
        data.use_curvature_filter = 'false';
        
        data = addOptimizerSettings(data, S.optRRA);
        data.use_verbose_printing = 'false';
             
    
        
    % COMPUTED MUSCLE CONTROL (CMC)
    % -----------------------------       
    case 'cmc',

        % CMCTool
        root = 'CMCTool';
        data.ATTRIBUTE.name = F.C3DFILE;
        
        data.model_file = F.scaledAdjustedModelFile;
        data.replace_actuator_set = 'false';
        data.actuator_set_files = F.actuatorSetFile.CMC;
        data.results_directory = F.CMC_Recults;
        data = addInitialFinalTimes(data, timeRange);
        
        data.compute_average_residuals = 'true';
        data.adjust_com_to_reduce_residuals = 'false';
        
        data.desired_kinematics_file = F.IKFile;
        data.lowpass_cutoff_frequency = -10;
        data = addExternalLoadSettings(data, S, F);
        data = addIntegratorSettings(data, S.inte);

        data.task_set_file = F.taskSetFile.CMC;
        data.constraints_file = F.controlConstraintsFile.CMC;       
        data.cmc_time_window = 0.01;
        data.use_curvature_filter = 'false';
        
        data = addOptimizerSettings(data, S.optCMC);
        data.use_verbose_printing = 'false';
        
        
        % CMCTool -> AnalysisSet
        data.AnalysisSet.ATTRIBUTE.name = 'Analyses';
        data = addKinematicsAnalysis(data);
        data = addActuationAnalysis(data);
        data = addBodyKinematicsAnalysis(data);



    
    % FORWARD DYNAMICS (FORWARD)
    % -----------------------------       
    case 'forward',
        
        % ForwardTool
        root = 'ForwardTool';
        data.ATTRIBUTE.name = F.C3DFILE;
        data.model_file = F.scaledAdjustedModelFile;
        data.replace_actuator_set = 'false';
        data.actuator_set_files = F.actuatorSetFile.CMC;
        data.results_directory = 'Forward_Results';
        data = addInitialFinalTimes(data, timeRange);
        data = addIntegratorSettings(data, S.inte);
        
        data.controls_file = F.controlsFile;
        data.states_file = F.statesFile;
        data.use_specified_dt = 'true';
        data = addExternalLoadSettings(data, S, F);
        
        % ForwardTool -> AnalysisSet
        data.AnalysisSet.ATTRIBUTE.name = 'Analyses';
        data = addKinematicsAnalysis(data);
        data = addActuationAnalysis(data);
        data = addBodyKinematicsAnalysis(data);
        
        
       
    % PETURB ANALYSIS (NOT ADDED YET)
    % --------------------------------         
           
        
    
        
    % UNKNOWN XML TYPE
    % ----------------
    otherwise
        fprintf('Invalid XML type. The type is case sensitve\n%s', types);
end



% Save the data structure to an xml file
% --------------------------------------
if saveFile > 0
    fileName = sprintf('%s_Setup_%s.xml', F.C3DFILE, upper(type));
    xmlString = xml_formatany(data, root);

    fid = fopen(fileName, 'w');
    if fid < 0
        fprintf('\nERROR: %s could not be opened for writing...\n\n', fileName);
        return
    end
    fprintf(fid, '%s\n', xmlString);
    fclose(fid);

    % xml_save(fileName, data, 'off');
    fprintf('%s saved sucessfully...\n', fileName);
end







% ========================================================================
% SUBFUNCTION: data = addIntegratorSettings(data, inte)
% ========================================================================
% Loads the data structure with integrator information
% Used for RRA1, RRA2, CMC and FORWARD

function data = addIntegratorSettings(data, inte)

data.maximum_number_of_integrator_steps = inte.maximum_number_of_integrator_steps;
data.maximum_integrator_step_size = inte.maximum_integrator_step_size;
data.integrator_error_tolerance = inte.integrator_error_tolerance;
data.integrator_fine_tolerance = inte.integrator_fine_tolerance;
data.output_precision = inte.output_precision;

  




% ========================================================================
% SUBFUNCTION: data = addExternalLoadSettings(data, S, F)
% ========================================================================
% Loads the data structure with "experimental" external load information
% Used for ID, RRA1, RRA2, CMC and FORWARD
% Note: If kinetics are modified in RRA for example, run this function
%       and then overwrite data.external_loads_file with the new source

function data = addExternalLoadSettings(data, S, F)

data.external_loads_file = F.IKFile;
data.external_loads_model_kinematics_file = F.IKFile;
data.external_loads_body1 = 'calcn_r';      % Right foot first
data.external_loads_body2 = 'calcn_l';      % Left foot second
data.lowpass_cutoff_frequency_for_load_kinematics = S.kinematics_filter_freq;






% ========================================================================
% SUBFUNCTION: data = addOptimizerSettings(data, opt)
% ========================================================================
% Loads the data structure with optimizer information
% Used for RRA1, RRA2 and CMC

function data = addOptimizerSettings(data, opt)

data.use_fast_optimization_target = opt.use_fast_optimization_target;
data.optimizer_derivative_dx = opt.optimizer_derivative_dx;
data.optimizer_convergence_criterion = opt.optimizer_convergence_criterion;
data.optimizer_max_iterations = opt.optimizer_max_iterations;
data.optimizer_print_level = opt.optimizer_print_level;
data.optimizer_algorithm = opt.optimizer_algorithm;






% ========================================================================
% SUBFUNCTION: data = addInitialFinalTimes(data, C3Dkey)
% ========================================================================
% Loads the data structure with initial and final time information
% Used for ID, RRA1, RRA2, CMC and FORWARD

function data = addInitialFinalTimes(data, timeRange)

data.initial_time = timeRange(1);
data.final_time = timeRange(end);






% ========================================================================
% SUBFUNCTION: data = addKinematicsAnalysis(data)
% ========================================================================
% Adds a KINEMATICS analysis

function data = addKinematicsAnalysis(data)

data.AnalysisSet.objects.Kinematics.ATTRIBUTE.name = 'Kinematics';
data.AnalysisSet.objects.Kinematics.on = 'true';
data.AnalysisSet.objects.Kinematics.step_interval = 10;
data.AnalysisSet.objects.Kinematics.in_degrees = 'true';






% ========================================================================
% SUBFUNCTION: data = addActuationAnalysis(data)
% ========================================================================
% Adds an ACTUATION analysis

function data = addActuationAnalysis(data)

data.AnalysisSet.objects.Actuation.ATTRIBUTE.name = 'Actuation';
data.AnalysisSet.objects.Actuation.on = 'true';
data.AnalysisSet.objects.Actuation.step_interval = 10;
data.AnalysisSet.objects.Actuation.in_degrees = 'true';






% ========================================================================
% SUBFUNCTION: data = addBodyKinematicsAnalysis(data)
% ========================================================================
% Adds a BODY KINEMATICS analysis

function data = addBodyKinematicsAnalysis(data)

data.AnalysisSet.objects.BodyKinematics.ATTRIBUTE.name = 'BodyKinematics';
data.AnalysisSet.objects.BodyKinematics.on = 'true';
data.AnalysisSet.objects.BodyKinematics.step_interval = 10;
data.AnalysisSet.objects.BodyKinematics.in_degrees = 'true';






% ========================================================================
% SUBFUNCTION: data = addInverseDynamicsAnalysis(data)
% ========================================================================
% Adds an INVERSE DYNAMICS analysis

function data = addInverseDynamicsAnalysis(data)

data.AnalysisSet.objects.InverseDynamics.ATTRIBUTE.name = 'InverseDynamics';
data.AnalysisSet.objects.InverseDynamics.on = 'true';
data.AnalysisSet.objects.InverseDynamics.step_interval = 1;
data.AnalysisSet.objects.InverseDynamics.use_model_actuator_set = 'false';

