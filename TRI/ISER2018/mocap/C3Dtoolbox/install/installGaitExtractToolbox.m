% Install script for GaitExtract Toolbox
% Tim Dorn
% 28th June 2008
% --------------------------------------------------------------
% 
% Instructions
% ------------
% 
% 1) Copy the ENTIRE toolbox directory to a desired location
% 2) Load MATLAB and change directory to this install directory
% 3) Run this install program
% 
% --------------------------------------------------------------

fprintf('Instructions\n')
fprintf('------------\n')
fprintf('1) Copy the ENTIRE toolbox directory to a desired location\n')
fprintf('2) Load MATLAB and change directory to this install directory\n')
fprintf('3) Run this install program\n\n')

cd ..
p = pwd;
cd install
fprintf('Toolbox will be added to the matlab path in: %s\n', p)

R = input('Ready to proceed with installation? (y/n): ', 's');
if R ~= 'y',
    fprintf('Installation Halted...\n')
    return
end


% Step 1: Install c3dServer (Motion Labs)
% ---------------------------------------
fprintf('Step 1: Install c3dServer & C3DViewer (Motion Labs)\n')
!c3dserver_install.exe
!mlsviewer_install.exe


% Step 2: Set up matlab paths
% ---------------------------
fprintf('Step 2: Set up Toolbox Paths\n')
pathLine = sprintf('addpath(genpath(''%s''))', p);
eval(pathLine)
savepath
fprintf('Installed Sucessfully!\nToolbox is ready to use RIGHT NOW!\n')

