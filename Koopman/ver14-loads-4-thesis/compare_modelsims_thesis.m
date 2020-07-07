% compare_modelsims_thesis
%   Simulates linear, bilinear, and nonlinear realizations of same system
%   and compares them to data from a validation data trial
%
%   This is not a general purpose function.

%% Choose whether to save siulation results or not
saveon = true;

%% Load the models

% load the Arm system class
temp_arm = load('C:\Users\danie\OneDrive\Documents\MATLAB\UM-Matlab\Koopman\ver14-loads-4-thesis\systems\thesis-arm-markers_noload_3-mods_1-links_20hz\thesis-arm-markers_noload_3-mods_1-links_20hz.mat');
Arm = temp_arm.Arm;

% load the model classes
temp_lin = load('C:\Users\danie\OneDrive\Documents\MATLAB\UM-Matlab\Koopman\ver14-loads-4-thesis\systems\thesis-arm-markers_noload_3-mods_1-links_20hz\models\linear_poly-3_n-6_m-3_del-0_2020-06-09_16-42.mat');
temp_bilin = load('C:\Users\danie\OneDrive\Documents\MATLAB\UM-Matlab\Koopman\ver14-loads-4-thesis\systems\thesis-arm-markers_noload_3-mods_1-links_20hz\models\bilinear_poly-3_n-6_m-3_del-0_2020-06-09_16-43.mat');
temp_nonlin = load('C:\Users\danie\OneDrive\Documents\MATLAB\UM-Matlab\Koopman\ver14-loads-4-thesis\systems\thesis-arm-markers_noload_3-mods_1-links_20hz\models\nonlinear_poly-3_n-6_m-3_del-0_2020-06-13_14-10.mat');
Ksysid_lin = temp_lin.sysid_class;  % linear model
Ksysid_bilin = temp_bilin.sysid_class;  % bilinear model
Ksysid_nonlin = temp_nonlin.sysid_class;  % nonlinear model

% load the real system trajectory
valres = load('C:\Users\danie\OneDrive\Documents\MATLAB\UM-Matlab\Koopman\ver14-loads-4-thesis\systems\thesis-arm-markers_noload_3-mods_1-links_20hz\simulations\tf-30s_ramp-06s.mat');

%% Simulate the response of each model

linres = Ksysid_lin.val_model( Ksysid_lin.model , Ksysid_lin.valdata{1} );
bilinres = Ksysid_bilin.val_BLmodel( Ksysid_bilin.model , Ksysid_lin.valdata{1} );
nonlinres = Ksysid_nonlin.val_NLmodel( Ksysid_nonlin.model , Ksysid_lin.valdata{1} );

%% Plot the results

figure; % linear model
subplot(6,1,1)
hold on;
plot( linres.real.t , linres.real.y(:,1) , 'b' );
plot( linres.sim.t , linres.sim.y(:,1) , 'r' );
hold off;
subplot(6,1,2)
hold on;
plot( linres.real.t , linres.real.y(:,2) , 'b' );
plot( linres.sim.t , linres.sim.y(:,2) , 'r' );
hold off;
subplot(6,1,3)
hold on;
plot( linres.real.t , linres.real.y(:,3) , 'b' );
plot( linres.sim.t , linres.sim.y(:,3) , 'r' );
hold off;
subplot(6,1,4)
hold on;
plot( linres.real.t , linres.real.y(:,4) , 'b' );
plot( linres.sim.t , linres.sim.y(:,4) , 'r' );
hold off;
subplot(6,1,5)
hold on;
plot( linres.real.t , linres.real.y(:,5) , 'b' );
plot( linres.sim.t , linres.sim.y(:,5) , 'r' );
hold off;
subplot(6,1,6)
hold on;
plot( linres.real.t , linres.real.y(:,6) , 'b' );
plot( linres.sim.t , linres.sim.y(:,6) , 'r' );
hold off;


figure; % bilinear model
subplot(6,1,1)
hold on;
plot( bilinres.real.t , bilinres.real.y(:,1) , 'b' );
plot( bilinres.sim.t , bilinres.sim.y(:,1) , 'r' );
hold off;
subplot(6,1,2)
hold on;
plot( bilinres.real.t , bilinres.real.y(:,2) , 'b' );
plot( bilinres.sim.t , bilinres.sim.y(:,2) , 'r' );
hold off;
subplot(6,1,3)
hold on;
plot( bilinres.real.t , bilinres.real.y(:,3) , 'b' );
plot( bilinres.sim.t , bilinres.sim.y(:,3) , 'r' );
hold off;
subplot(6,1,4)
hold on;
plot( bilinres.real.t , bilinres.real.y(:,4) , 'b' );
plot( bilinres.sim.t , bilinres.sim.y(:,4) , 'r' );
hold off;
subplot(6,1,5)
hold on;
plot( bilinres.real.t , bilinres.real.y(:,5) , 'b' );
plot( bilinres.sim.t , bilinres.sim.y(:,5) , 'r' );
hold off;
subplot(6,1,6)
hold on;
plot( bilinres.real.t , bilinres.real.y(:,6) , 'b' );
plot( bilinres.sim.t , bilinres.sim.y(:,6) , 'r' );
hold off;


figure; % bilinear model
subplot(6,1,1)
hold on;
plot( nonlinres.real.t , nonlinres.real.y(:,1) , 'b' );
plot( nonlinres.sim.t , nonlinres.sim.y(:,1) , 'r' );
hold off;
subplot(6,1,2)
hold on;
plot( nonlinres.real.t , nonlinres.real.y(:,2) , 'b' );
plot( nonlinres.sim.t , nonlinres.sim.y(:,2) , 'r' );
hold off;
subplot(6,1,3)
hold on;
plot( nonlinres.real.t , nonlinres.real.y(:,3) , 'b' );
plot( nonlinres.sim.t , nonlinres.sim.y(:,3) , 'r' );
hold off;
subplot(6,1,4)
hold on;
plot( nonlinres.real.t , nonlinres.real.y(:,4) , 'b' );
plot( nonlinres.sim.t , nonlinres.sim.y(:,4) , 'r' );
hold off;
subplot(6,1,5)
hold on;
plot( nonlinres.real.t , nonlinres.real.y(:,5) , 'b' );
plot( nonlinres.sim.t , nonlinres.sim.y(:,5) , 'r' );
hold off;
subplot(6,1,6)
hold on;
plot( nonlinres.real.t , nonlinres.real.y(:,6) , 'b' );
plot( nonlinres.sim.t , nonlinres.sim.y(:,6) , 'r' );
hold off;
