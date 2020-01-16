function ref_spiral = spiral_generation(starting_point,revolution,wide,height)
%staring_point is the intial x,y,z cooridnates of the spiral, where the
%direction of each axis is consistent with the one in mocap. (Unit:mm)

%revolution is an integer input that specifies the revolution of the
%whole spiral

%wide and height are the two coverage limits for spiral. Unit(:mm)

%ref_spiral is an N by 3 array.

    
      
    % starting_point = [0,0,0];
    % revolution = 2;
    % wide = 200;
    % height = 80;

    t = 0 : pi/20 : revolution*2*pi;
    a = linspace(0,wide/2,length(t));
    h = linspace(0,height,length(t));

    x = (a.*sin(t))';
    z = (a.*cos(t))';
    y = h';

    ref_spiral = [x,y,z] + starting_point;

    % figure;
    % plot(x,z)
    % axis equal
    % 
    % figure;
    % plot3(ref(:,1),ref(:,3),ref(:,2));

end

