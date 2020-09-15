% compare_modelpredictions_thesis
%   Simulates linear, bilinear, and nonlinear realizations of same system
%   and compares them to data from a validation data trial
%
%   This is not a general purpose function.

%% Choose whether to save siulation results or not
saveon = true;

%% Load the models

% load the Arm system class
temp_arm = load('C:\Users\danie\OneDrive\Documents\MATLAB\UM-Matlab\Koopman\ver14-loads-4-thesis\systems\thesis-arm-markers_noload_3-mods_1-links_20hz\thesis-arm-markers_gravload1_3-mods_1-links_20hz.mat');
Arm = temp_arm.Arm;

% load the model classes
temp_lin = load('C:\Users\danie\OneDrive\Documents\MATLAB\UM-Matlab\Koopman\ver14-loads-4-thesis\systems\thesis-arm-markers_noload_3-mods_1-links_20hz\models\linear_poly-3_n-6_m-3_del-0_2020-06-09_16-42.mat');
temp_bilin = load('C:\Users\danie\OneDrive\Documents\MATLAB\UM-Matlab\Koopman\ver14-loads-4-thesis\systems\thesis-arm-markers_noload_3-mods_1-links_20hz\models\bilinear_poly-3_n-6_m-3_del-0_2020-06-09_16-43.mat');
temp_nonlin = load('C:\Users\danie\OneDrive\Documents\MATLAB\UM-Matlab\Koopman\ver14-loads-4-thesis\systems\thesis-arm-markers_noload_3-mods_1-links_20hz\models\nonlinear_poly-3_n-6_m-3_del-0_2020-06-13_14-10.mat');
Ksysid_lin = temp_lin.sysid_class;  % linear model
Ksysid_bilin = temp_bilin.sysid_class;  % bilinear model
Ksysid_nonlin = temp_nonlin.sysid_class;  % nonlinear model
