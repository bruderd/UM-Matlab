% sim_debug.m
% 
% Quick and dirty simulation file for debugging the system.
% Run setup.m before executing this script, sys is needed.

params = sys.params;

%% set constant input and initial conditions

u = pi/6 * [1 0.5 -1]' ; %3 * ones( params.Nmods , 1);

a0 = zeros( params.Nlinks , 1 );
adot0 = zeros( params.Nlinks , 1 );

%% simulate system
tic;    % for measureing elapsed time

tf = 1;    % final time
nsteps = 3000;
tspan = linspace( 0 , tf , nsteps * tf );
% [ t , y ] = ode45( @(t,x) vf_sym(t,x,u) , tspan , [ a0 ; adot0 ] );  % with symbolic inversion
% [ t , y ] = ode45( @(t,x) vf_sys(t,x,u,sys) , tspan , [ a0 ; adot0 ] );  % with numerical inversion
[ y ] = ode5( @(t,x) vf_sys(t,x,u,sys) , tspan , [ a0 ; adot0 ] );  % with numerical inversion, fixed time step
t = tspan';

simTime = toc;  % measure time for the simulation
%% plot steady state robot configuration

alpha = y(: , 1:params.Nlinks );   % joint angles over time

% plot alpha over time
figure; plot(t,alpha);

% plot final robot configuration
figure;
plot_manip( alpha(end,:)' , params );

%% animate the result

fig = figure;   % create figure for the animation
axis([-params.L, params.L, -0.5*params.L, 1.5*params.L])
set(gca,'Ydir','reverse')
xlabel('x(m)')
ylabel('y(m)')

% Prepare the new file.
vidObj = VideoWriter( ['animations' , filesep , params.sysName , '.mp4'] , 'MPEG-4' );
vidObj.FrameRate = 30;
open(vidObj);

set(gca,'nextplot','replacechildren', 'FontUnits' , 'normalized');

totTime = tf;    % total time for animation (s)
totFrames = 30 * totTime;   % total frames in 30 fps video
k = 1;  % initialize counter



for i = 1:totFrames
    
    index = (i-1) * floor( nsteps / totFrames ) + 1;   % skips points between frames
    
    [ X , ~ ] = alpha2x( alpha(index,:)', params );
    x = [0; X(1:2:end)];
    y = [0; X(2:2:end)];
    marker = sys.get_markers( alpha(index,:) );   % get mocap sensor location
    [shape , ~ ] = sys.get_shape( alpha(index,:) , 5); % get polynomial approx of shape (3rd order)
    
    hold on;
    p1 = plot(x, y, 'b-o');
    p2 = plot( marker(:,1) , marker(:,2) , 'r*');
    p3 = plot( shape(:,1) , shape(:,2) , 'r');
    hold off;
    
    % write each frame to the file
    currFrame = getframe(fig);
    writeVideo(vidObj,currFrame);
    
    delete(p1); delete(p2); delete(p3);
    
end

close(vidObj);
