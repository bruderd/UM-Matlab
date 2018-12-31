% manipulator_main.m
%   Main script for simulating a double pendulum system

%% set manipulator parameters
params = struct;

% system parameters
params.name = 'manipulator_1000s';
params.ploton = true;       % decides whether to animate the result
params.movie = false;       % decides wheter to save animation or not
params.n = 4;   % number of states
params.p = 2;   % number of inputs
params.ny = 2;  % dimension of the output

% physical parameters (note, angles are absolute, not relative)
params.phi1                = 0; % (pi/2)*rand - pi/4;
params.dtphi1              = 0;
params.phi2                = 0; % (pi/2)*rand - pi/4;
params.dtphi2              = 0;
params.g                   = 9.81; 
params.m1                  = 1; 
params.m2                  = 1; 
params.l1                  = 1; 
params.l2                  = 1;
params.k1                  = 50;
params.k2                  = 50;
params.b1                  = -3.5;
params.b2                  = -3.5;

% simulation parameters
params.Ts = 0.02;
params.tf = 1000;
params.x0 = [params.phi1, params.phi2, params.dtphi1, params.dtphi2]';

% step/ramp input parameters
params.steplen = 1;  % duration of each step
params.steps = 2*pi * ( rand(2000, params.p) - 0.5 );   % random list of 1000 step inputs between -pi and pi
params.tconstant = 1.5;

%% simulate manipulator
data = gen_fakeData( params.name, @manipulator_vf, @manipulator_input, params );

%% calculate output, the position of the end effector

data.y = zeros(size(data.x));
data.y(:,1) = params.l1 * sin( data.x(:,1) ) + params.l2 * sin( data.x(:,2) );
data.y(:,2) = params.l1 * cos( data.x(:,1) ) + params.l2 * cos( data.x(:,2) );

%% save function that evaluates dynamics with these system parameters
x = sym('x' , [params.n , 1]);
u = sym('u' , [params.p , 1]);
symbolic_dynamics = manipulator_vf(x, u, params);


%% show animation of the system
fig = figure;
h=plot(0,0,'MarkerSize',30,'Marker','.','LineWidth',2);
range=1.1*(params.l1+params.l2); axis([-range range -range range]); axis square;
set(gca,'nextplot','replacechildren');
    for i=1:length(data.x(:,1))-1
        if (ishandle(h)==1)
            Xcoord=[0, params.l1*sin(data.x(i,1)) , params.l1*sin(data.x(i,1))+params.l2*sin(data.x(i,2)) ];
            Ycoord=[0 ,-params.l1*cos(data.x(i,1)) , -params.l1*cos(data.x(i,1))-params.l2*cos(data.x(i,2)) ];
            set(h,'XData',Xcoord,'YData',Ycoord);
            title( [ 't = ' , num2str(data.t(i)) ] );
            drawnow;
            F(i) = getframe;
            if params.movie==false
                pause(data.t(i+1)-data.t(i));
            end
        end
    end
    if params.movie==true
        movie2avi(F,'doublePendulumAnimation.avi','compression','Cinepak','fps',fps)
    end
    

%% save y as x so that it can work with the sysid code
data.



