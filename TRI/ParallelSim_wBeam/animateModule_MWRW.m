function animateModule_MWRW( t, y, u, params )
% animateModule.m
%   Animates the end effector of a module of parallel FREEs driven with
%   sinusoidal pressure inputs slightly out of phase with one another.

% dt = 0.01;   % size of time step
% T = 2*pi;     % final time
% t = 0:dt:T; % time vector
% timetodelay = dt;

filename = 'MWRWsim0.gif';

fig = figure('Position', [1 41 1536 748.8000]); % make the figure bigger

% interpolate so that time(s) where figure is drawn are equally spaced
tq = linspace(0,t(end),params.frames)';
yq = interp1(t,y,tq);
uq = interp1(t,u,tq);

for i = 1:length(tq)
    if i == 1
        timetodelay = tq(1) * 1e0;
    else
        timetodelay = (tq(i)-tq(i-1)) * 1e0;
    end

    % do animation
    clf
    drawModule_MWRW(yq(i,:)', params, tq(i), uq(i,:));
    
    % create a gif of animation
    frame = getframe(fig);  % take a snapshot of the figure
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1
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