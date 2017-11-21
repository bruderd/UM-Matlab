function params = setParams(p, n , module, free, sim)
%setParams: Creates parameter struct which holds all parameters related to
%soft manipulator.
%   Inputs:
%       p: number of modules in manipulator
%       n: number of actuators in each module, must be vector of length p
%       module: all parameters related to each module
%       free: all paramters that are related to each actuator
%       sim: simulation parameters

params = struct;

%% Make sure the sizes of inputs are consistent
if size(n,1) ~= p
    error('n must have p rows');
elseif size(module,1) ~= p
    error('Not enough/too many module parameters');
elseif size(free,1) ~= sum(n)
    error('Not enough/too many free parameters');
end


%% MANIPULATOR PARAMETERS
params.p = p;
params.n = n;

%% MODULE PARAMETERS
% module = [L, block density, EI of spine], size(module) = p x 3

params.L = module(:, 1);
params.denBlock = module(:, 2);
params.dimBlock = module(:, 3:5);
params.EIspine = module(:, 6);   

% assuming each module connection block is a rectangular prism
params.m = params.denBlock .* params.dimBlock(:,1) .* params.dimBlock(:,2) .* params.dimBlock(:,3);     % mass of the module end blocks (kg)
params.I = zeros(3*p, 3);   % moment of inertia matrix blocks in local coordinate frame
for i = 1:p
    params.I(3*(i-1)+1 : 3*i, :) = [params.m(i)/12 * (params.dimBlock(i,2)^2 + params.dimBlock(i,3)^2), 0, 0;...
                                    0, params.m(i)/12 * (params.dimBlock(i,1)^2 + params.dimBlock(i,3)^2), 0;...
                                    0, 0, params.m(i)/12 * (params.dimBlock(i,1)^2 + params.dimBlock(i,2)^2)];
end


%% FREE PARAMETERS
% free = [Gama, R, xattach, yattach], size(free) = sum(n) x 4

params.Gama = free(:,1);
params.R = free(:,2);
params.attach = free(:,3:4);
params.xattach = free(:,3);
params.yattach = free(:,4);

% the length of each FREE (determined by length of module it's in)
params.Lfree(1:n(1),1) = params.L(1);
for i = 2:p
    params.Lfree(sum(n(1:i-1))+1 : sum(n(1:i)), 1) = params.L(i);
end

params.B = abs(params.Lfree ./ cos(params.Gama));   % fiber length (must be positive))
params.Nf = -params.Lfree ./ (2*pi*params.R) .* tan(params.Gama); % total fiber windings in revolutions (when relaxed)


%% SIMULATION PARAMETERS
% EDIT THIS LATER WHEN YOU ARE ACTUALLY GOING TO RUN SIMULATION

% initial conditions of simulation
% params.X0_0 = zeros(2*6*p, 1);
% for i = 1:p
%     params.X0_0(2*6*(i-1)+3, 1) = sum(params.L(1:i));
% end
params.X0dot_0 = zeros(2*3*p, 1);
params.X0_0 = zeros(2*3*p, 1);

% % Initial conditions (could be made more generic in the future)
% params.x0 = [0, 0, 0, 0, 0, 0]';
% params.xdot0 = [params.x0(4), params.x0(5), params.x0(6), 0, 0, 0]';
% params.u0 = [1000, 0, 1000, 0]';
% 
% params.T = 1;   %final time
% params.N = 50; %number of steps
% params.dt = params.T/params.N;    %size of one time step
% params.n = length(params.x0);  %dimension of state vector x
% params.m = length(params.u0);  %dimension of state vector u
% 
% % Maximum pressure for the FREEs
% params.Pmax = 400e3;    % (Pa)


%% VISUALISATION PARAMETERS





        


end

