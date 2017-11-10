% RotaryPendulum_MLdemo    This demo illustrates an application of the
%                Euler-Lagrange equation. The system of differential
%                equations derived with the EulerLagrange tool is solved
%                using ode45. 
%                Note: If the Simulink 3D Animation Toolbox is available 
%                the dynamics of the system can be visualized in the VR
%                world.
%
% Copyright 2015-2016 The MathWorks, Inc.

% Solve DE numerically using ode45
% Set simulation time and initial state vector
tspan = [0 10];
X0    = [0 0 10 0]*pi/180;
% Set parameter values
I     = 1;
m     = 1;
l     = 1;
R     = 1;
g     = 9.81;
tau   = 0;    % No external input torque

% Use created .m file to solve DE
[t, X]  = ode45(@RotaryPendulum_ODE,tspan,X0,[],tau,I,m,l,R,g);

% Plot state vector
plot(t,X)
title('Rotary pendulum')
xlabel('t')
ylabel('X(t)')
legend('{\alpha}','d{\alpha}/dt','{\beta}','d{\beta}/dt')

% Start animation if VR toolbox is available
hasVR = license('test', 'virtual_reality_toolbox');
if ~hasVR
    msgbox('You do not seem to have the Simulink 3D Animation Toolbox.',...
           'Toolbox missing');
else
    % Open the VirtualReality world and display in figure
    world = vrworld('RotaryPendulum_VR.wrl');
    open(world);
    fig   = view(world, '-internal');
    
    % Define VR nodes to be able to rotate objects
    RotaryArm = vrnode(world, 'RotaryArm');
    Pendulum  = vrnode(world, 'Pendulum');
    
    % Loop through state vector and re-draw VR
    for ii = 1:size(X,1)
        RotaryArm.rotation = [0 1 0 X(ii,1)];
        Pendulum.rotation  = [0 0 1 X(ii,3)];
        vrdrawnow;
        pause(0.05);
    end
end % if statement
% End of script