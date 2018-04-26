function params = setParams(filename)
%setParams: Creates a struct containing all parameters of the problem
%   User defines parameter values in the first section, then all parameters
%   are stored in the struct 'params'.

%% (USER EDIT) User defined parameters

% Actuator parameters
num = 3;    % number of FREEs in combination
Gama = deg2rad([15, -15, 52.77]); % relaxed fiber angle of each FREE
R = (12.319e-3)/2 * ones(1,num);  % relaxed radius of each FREE [m]
L = 0.16 * ones(1,num);   %  relaxed length of each FREE [m] 
d = [0, -0.02286, 0 ; 0.01979676, 0.01143, 0 ; -0.01979676, 0.01143, 0]'; % location of attachment points to the end effector [m]
a = [0,0,1 ; 0,0,1 ; 0,0,1]';    % direction of FREE axis at attachment point [unit vector]
pmin = (1/0.145038) * 1e3 * [1 1 1];   % min gauge pressure for each FREE [Pa]
pmax = (1/0.145038) * 1e3 * [12 12 20];   % max gauge pressure for each FREE [Pa]

% End effector parameters
deff = [0,0,0]; % location of origin of end effector coordinates in global coordinates
euleff = eye(3);  % orientation of end effector coordinate frame in global coordinates, written as rotation matrix
meff = 0;   % mass of the end effector [kg]
cmeff = [0,0,0]';   % location of the center of mass of end effector [m]
C = -(1)*[1e1 0 0 1e-3; 1e1 0 0 1e-3; 1e1 0 0 1e-3]';   % compliance (stiffness) matrix for each FREE vectorized so that [c1, c2; c3, c4] = [c1, c2, c3, c4]', horizontally concatenated

% QP parameters
tol = 1e-1;   % constraint tolerance of the QP

% SysID parameters
psteps = 8;     % how finely to break up pmax

% Enfield TR parameters
TRpsimax = 25;      % pressure (in psi) that corresponds to 10V input signal to TR pressure regulator

% Mocap parameters (may or may not use...)
topIDs = [9, 13, 15];    % marker id numbers for LEDs on the top block
effIDs = [8, 10, 12, 14];  % makrer id numbers for LEDs on the bottom block

%% check that the sizes of parameters entered are consistent
if ~(all(size(L) == size(R)) && all(size(R) == size(Gama)) && all(size(Gama) == size(pmin))...
        && all(size(pmin) == size(pmax)) && all(size(pmax) == size(d(1,:))) && all(size(d) == size(a)))
    error('The sizes of one or more assigned variables are not consistent');
end
    
%% Create struct to store all parameters
params = struct;

params.num = num;
params.Gama = Gama;
params.R = R;
params.L = L;
params.d = d;
params.a = a;
params.pmin = pmin;
params.pmax = pmax;
params.meff = meff; 
params.cmeff = cmeff;   
params.deff = deff;
params.euleff = euleff;
params.psteps = psteps;
params.TRpsimax = TRpsimax;
params.topIDs = topIDs;
params.effIDs = effIDs;
params.tol = tol;

params.B = abs(params.L ./ cos(params.Gama));   % fiber length (must be positive))
params.N = -params.L ./ (2*pi*params.R) .* tan(params.Gama); % total fiber windings in revolutions (when relaxed)

% force transformation matrix from FREE to end effector coordinates
% (cumulative)
params.D = zeros(6,2*num);
for i = 1:num
    dix = [0, -d(3,i), d(2,i); d(3,i), 0, -d(1,i); -d(2,i), d(1,i), 0];
    params.D(: , 2*(i-1)+1:2*i) = [[a(:,i), zeros(3,1)] ; [flipud(dix*a(:,i)), zeros(3,1)] + [zeros(3,1), flipud(a(:,1))]];  % Di's are horizontally concatenated
end

% compliance matrix (cumulative)
params.C = zeros(2*num, 2*num);
for i = 1:num
    params.C(2*(i-1)+1 : 2*i, 2*(i-1)+1 : 2*i) = reshape( C(:,i), [2,2] )';
end

% penalty weighting (this is used to focus on the equilibrium point with lowest pressure)
params.penalty = 1e-5;
params.tol = 1e-2;  % tolerance to within which equality constraint of quadprog will  be satisfied

%% set the inverse kinematic relationship for the system
params = setInvKin(params); % a few parameters are added in this function

%% save these parameters as a .mat file

% check for optional argument, if given, save params as .mat file with that name
if exist('filename','var')
    current_folder = cd;
    savetolocation = strcat(current_folder, '\configs\', filename);
    save(char(savetolocation), 'params');
end

end




