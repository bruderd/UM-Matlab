function animate_doublePendulum( sol_real, sol_sysid, params )
%double_pendulum_animate: animates double pendulum simulation
%   sol :       output of ode45 solution
%   params :    parameters of the system

% Extract inputs from params struct
duration = params.duration;
fps = params.fps;
movie = params.movie;

nframes=duration*fps;
t = linspace(0,duration,nframes);
yr = deval(sol_real,t);
% ys = deval(sol_sysid,t);
ys = interp1(sol_sysid.x, sol_sysid.y', t)';

phi1r=yr(1,:)'; % dtphi1r=yr(3,:)';
phi2r=yr(2,:)'; % dtphi2r=yr(4,:)';
phi1s=real(ys(1,:)'); % dtphi1s=ys(3,:)';
phi2s=real(ys(2,:)'); % dtphi2s=ys(4,:)';
l1=params.l1; l2=params.l2;

figure;
hold on
hr=plot(0,0,'MarkerSize',30,'Marker','.','LineWidth',2);
hs=plot(0,0,'MarkerSize',30,'Marker','.','LineWidth',2, 'color', 'r');
range=1.1*(l1+l2); axis([-range range -range range]); axis square;
set(gca,'nextplot','replacechildren');

for i=1:length(phi1r)-1
    if (ishandle(hs)==1)
        
        Xcoordr=[0,l1*sin(phi1r(i)),l1*sin(phi1r(i))+l2*sin(phi2r(i))];
        Ycoordr=[0,-l1*cos(phi1r(i)),-l1*cos(phi1r(i))-l2*cos(phi2r(i))];
        set(hr,'XData',Xcoordr,'YData',Ycoordr);
        
        % Coordinates if state is in terms of theta
        Xcoords=[0,l1*sin(phi1s(i)),l1*sin(phi1s(i))+l2*sin(phi2s(i))];
        Ycoords=[0,-l1*cos(phi1s(i)),-l1*cos(phi1s(i))-l2*cos(phi2s(i))];

        % Coordinates if state is in terms of xy
%         Xcoords=[0, ys(1,i), ys(3,i)];
%         Ycoords=[0, ys(2,i), ys(4,i)];
        
        % Coordinates if state is in terms of x2y2
%         Xcoords=[0, 0, ys(1,i)];
%         Ycoords=[0, 0, ys(2,i)];
        set(hs,'XData',Xcoords,'YData',Ycoords);
        
        drawnow;
        F(i) = getframe;
        if movie==false
            pause(t(i+1)-t(i));
        end
    end
end

hold off;

if movie==true
    vidObj = VideoWriter('doublePendulumAnimation', 'MPEG-4');
    vidObj.FrameRate = fps;
    open(vidObj);
    writeVideo(vidObj, F);
    close(vidObj);
end

end

