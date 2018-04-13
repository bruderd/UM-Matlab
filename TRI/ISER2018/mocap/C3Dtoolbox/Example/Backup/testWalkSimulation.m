% Example: Normal Level Walking Gait
% Muscle Actuated Analysis in OpenSim
% Tim Dorn
% -----------------------------------------------------------------------
% Note: This file can act as a template for extracting any experimental
% data from a motion capture system. loadlabels.m needs to be configured
% for the specific laboratory (see manual).
% -----------------------------------------------------------------------


%% Extract Data From Raw C3D Files
% ================================
clear all
close all

fprintf('\nAnswer YES when asked about the static trial.\n');


% File Name Descriptors
% ---------------------
c3dFileStatic = 'testStatic.c3d';
c3dFileDynamic = 'testWalking.c3d';


% Extract Event Keys
% ------------------
keySta = getEvents(c3dFileStatic, 2);
keyDyn = getEvents(c3dFileDynamic, 2, [1 2 3]);


% Extract Markers
% ---------------
markersSta = getMarkers(keySta, 'markersStatic');
markersDyn = getMarkers(keyDyn, 'markersDynamic');


% Extract Kinetics
% ----------------
opt = 2;        % Option to display all kinetic plots
[GRF, CoP, GRMo, GRMx] = getKinetics(keyDyn, opt, 0, markersDyn);


% Extract EMG
% -----------
[eVecGlob, EMGVecGlob] = batchEMGprocess(keyDyn, 'testEMG', 'EMGprocessTKE1', 'myEMG');
fprintf('\nLets plot some EMG signals...\n');
extractMotFile('file', 'testWalking_EMGset.mot');


% Write XML Files for Opensim
% ---------------------------
writeXML('scale', keySta);
writeXML('ik', keyDyn);


% Due to sensitivity to experimental data, the ID, RRA, and CMC XML
% setting files may need to be refined after they are written. The
% commented functions below output the general format. Parameters can then
% be changed (i.e. simulation times, integration parameters, etc). 
% The refined ID, RRA, and CMC files for this walking test case already
% exist in the example directory so make sure they are not overwritten.
% -----------------------------------------------------------------------
% writeXML('id', keyDyn);
% writeXML('rra', keyDyn);
% writeXML('cmc', keyDyn);


save extractedData

