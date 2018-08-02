function double_pendulum_animate( sol, params )
%double_pendulum_animate: animates double pendulum simulation
%   sol :       output of ode45 solution
%   params :    parameters of the system

% Extract inputs from params struct
duration = params.duration;
fps = params.fps;
movie = params.movie;

nframes=duration*fps;
t = linspace(0,duration,nframes);
y = deval(sol,t);

phi1=y(1,:)'; dtphi1=y(3,:)';
phi2=y(2,:)'; dtphi2=y(4,:)';
l1=params.l1; l2=params.l2;

figure;
h=plot(0,0,'MarkerSize',30,'Marker','.','LineWidth',2);
range=1.1*(l1+l2); axis([-range range -range range]); axis square;
set(gca,'nextplot','replacechildren');

    for i=1:length(phi1)-1
        if (ishandle(h)==1)
            Xcoord=[0,l1*sin(phi1(i)),l1*sin(phi1(i))+l2*sin(phi2(i))];
            Ycoord=[0,-l1*cos(phi1(i)),-l1*cos(phi1(i))-l2*cos(phi2(i))];
            set(h,'XData',Xcoord,'YData',Ycoord);
            drawnow;
            F(i) = getframe;
            if movie==false
                pause(t(i+1)-t(i));
            end
        end
    end
    if movie==true
        vidObj = VideoWriter('doublePendulumAnimation', 'MPEG-4');
        vidObj.FrameRate = fps;
        open(vidObj);
        writeVideo(vidObj, F);
        close(vidObj);
    end

end

