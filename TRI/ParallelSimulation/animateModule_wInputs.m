function animateModule_wInputs( t, y, u, params )
% animateModule.m
%   Animates the end effector of a module of parallel FREEs driven with
%   sinusoidal pressure inputs slightly out of phase with one another.

% dt = 0.01;   % size of time step
% T = 2*pi;     % final time
% t = 0:dt:T; % time vector
% timetodelay = dt;

filename = 'ModuleSimulation_wInputs2.gif';

fig = figure;

for i = 1:length(t)
    if i == 1
        timetodelay = t(1);
    else
        timetodelay = t(i)-t(i-1);
    end

    % do animation
    clf
    drawModule_wInputs(y(i,:)', params, t(i), u(i,:));
    
    % create a gif of animation
    frame = getframe(fig);  % take a snapshot of the figure
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1;
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf ,'DelayTime', timetodelay);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime', timetodelay);
    end
    
end





%% Check the driving pressure function
% figure
% hold on
% plot(t,P1)
% plot(t,P2)
% plot(t,P3)
% plot(t,P4)
% hold off