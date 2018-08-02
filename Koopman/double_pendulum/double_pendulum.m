function [sol, u] = double_pendulum(params)
% DOUBLE_PENDULUM Animates the double pendulum's (mostly) chaotic behavior.
%
%   author:  Alexander Erlich (alexander.erlich@gmail.com)
%
%   parameters:
%   
%   ivp=[phi1; dtphi1; phi2; dtphi2; g; m1; m2; l1; l2]
%
%                               Initial value problem. phi1 and dtphi1 are
%                               the initial angle and anglular velocity. g
%                               is gravity, m1 and l1 mass and rod length.
%                               For an explaining picture, see
%                               documentation file in same folder.
%  
%   duration                    The time interval on which the ode is
%                               solved spans from 0 to duration (in sec).
%
%   fps                         Frames Per Second. The framerate is
%                               relevant both for normal (realtime)
%                               animation and movie recording.
%
%   movie                       If false, a normal realtime animation of
%                               the motion of the double pendulum (the 
%                               framerate being fps) is shown.
%                               If true, a movie (.avi) is recorded. The
%                               filename is 'doublePendulumAnimation.avi'
%                               and the folder into which it is saved is
%                               the current working directory.
%
%   This function calls double_pendulum_ODE and is, in turn, called by
%   double_pendulum_init.
%
%   Example call:    >> double_pendulum([pi;0;pi;5;9.81;1;1;2;1],100,10,false)
%   Or, simply call  >> double_pendulum_init
%
%   ---------------------------------------------------------------------

clear All; clf;

% Extract inputs from params struct
ivp = [params.phi1; params.phi2; params.dtphi1; params.dtphi2; params.g;...
    params.m1; params.m2; params.l1; params.l2];
duration = params.duration;
fps = params.fps;
movie = params.movie;

% Simulate dynamics
sol = ode45(@(t, x) double_pendulum_ODE(t, x, get_u(t), params),[0 duration], ivp(1:4));
u = get_u(sol.x);       % get inputs corresponding to output points

% animate the results if desired
if movie==true
    double_pendulum_animate(sol, params);
end

%%  Old animation code 

% nframes=duration*fps;
% t = linspace(0,duration,nframes);
% y = deval(sol,t);
% u = get_u(t);
% 
% phi1=y(1,:)'; dtphi1=y(2,:)';
% phi2=y(3,:)'; dtphi2=y(4,:)';
% l1=ivp(8); l2=ivp(9);
% % phi1=x(:,1); dtphi1=x(:,2);
% % phi2=x(:,3); dtphi2=x(:,4);
% % l1=ivp(8); l2=ivp(9);
% 
% h=plot(0,0,'MarkerSize',30,'Marker','.','LineWidth',2);
% range=1.1*(l1+l2); axis([-range range -range range]); axis square;
% set(gca,'nextplot','replacechildren');
% 
%     for i=1:length(phi1)-1
%         if (ishandle(h)==1)
%             Xcoord=[0,l1*sin(phi1(i)),l1*sin(phi1(i))+l2*sin(phi2(i))];
%             Ycoord=[0,-l1*cos(phi1(i)),-l1*cos(phi1(i))-l2*cos(phi2(i))];
%             set(h,'XData',Xcoord,'YData',Ycoord);
%             drawnow;
%             F(i) = getframe;
%             if movie==false
%                 pause(t(i+1)-t(i));
%             end
%         end
%     end
%     if movie==true
%         vidObj = VideoWriter('doublePendulumAnimation', 'MPEG-4');
%         vidObj.FrameRate = fps;
%         open(vidObj);
%         writeVideo(vidObj, F);
%         close(vidObj);
%     end
    
    
    
    
    
    