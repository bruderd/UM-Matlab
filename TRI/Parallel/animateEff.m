% animateEff.m
%   Animates the end effector of a module of parallel FREEs driven with
%   sinusoidal pressure inputs slightly out of phase with one another.

dt = 0.1;   % size of time step
T = 2*pi;     % final time
t = 0:dt:T; % time vector

maxP = 1e6; % maximum FREE pressure

timetodelay = dt;
filename = 'parallelFREEs.gif';

fig = figure;

for i = 1:length(t)
%     % Pressures increase in magnitute over time
%     P1(i) = -cos(t(i))*maxP * log(t(i)+1)/2 + maxP;
%     P2(i) = -cos(t(i) + pi/2)*maxP * log(t(i)+1)/2 + maxP;
%     P3(i) = -cos(t(i) + 2*pi/2)*maxP * log(t(i)+1)/2 + maxP;
%     P4(i) = -cos(t(i) + 3*pi/2)*maxP * log(t(i)+1)/2 + maxP;

    % Pressures cycle and return to starting values
    P1(i) = -cos(t(i))*maxP + maxP;
    P2(i) = -cos(t(i) + pi/2)*maxP + maxP;
    P3(i) = -cos(t(i) + 2*pi/2)*maxP + maxP;
    P4(i) = -cos(t(i) + 3*pi/2)*maxP + maxP;
    
    P(:,i) = [P1(i); P2(i); P3(i); P4(i)];
  
    % do animation
    clf
    plotPos_iter(P(:,i), params);
    
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
figure
hold on
plot(t,P1)
plot(t,P2)
plot(t,P3)
plot(t,P4)
hold off