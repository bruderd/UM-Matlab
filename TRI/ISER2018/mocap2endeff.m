function [ endeff ] = mocap2endeff( mocap )
%mocap2x: Uses mocap marker position data to compute the state of the end
%   Detailed explanation goes here

%% calculate coordinates of topBlock and end effector (in mocap coordinates)
% calculate center of top block as the average position of top block markers over whole trial
total = zeros(1,3);
for i = 1 : size(mocap.topxyz, 3)
   total = total + nanmean(mocap.topxyz(:,:,i));
   topcent = total / i;
end

% calculate center of end effector as the average position of end effector markers at first 60 points (1 s)
total = zeros(1,3);
for i = 1 : size(mocap.topxyz, 3)
   total = total + nanmean(mocap.effxyz(1:60,:,i));
   effcent = total / i;
end

%% build transformation between coordinate systems
% end effector origin and coordinate axes expressed in mocap coordinates
origine = effcent;
ze = (effcent - topcent) / norm(effcent - topcent);
lede = nanmean(mocap.effxyz(1:60,:,1));    % pick one of the LEDs, i.e. the "first" one
xedir = (lede - effcent); % vector pointing from center to one of the LEDs, will point roughly along x-axis
ye = cross(ze,xedir) / norm( cross(ze,xedir) );
xe = cross(ye,ze) / norm( cross(ye,ze) );

% Homogeneous transformation from end effector coordinates to mocap coordinates
Re2m = [xe', ye', ze'];
He2m = [ Re2m, origine' ; [0 0 0 1] ];

% Homogeneous transformation from mocap coordinates to end effector coordinates
Rm2e = Re2m';
Hm2e = [ Rm2e, -Rm2e*origine' ; [0 0 0 1] ];

%% transform mocap data of end effector marker positions into end effector coordinates
endeff = struct;    % create output struct

for i = 1 : size(mocap.effxyz, 3)  % do for each LED mounted to end effector
   for j = 1 : mocap.frames
        jfoo = Hm2e*[mocap.effxyz(j,:,i),  1]';    % transform coordinates
        endeff.LEDxyz(j,:,i) = jfoo(1:3)';     % save in output struct
   end
end

%% calculate the state vector for the end effector at each frame

% postition of end effector, as average location of leds 
total = zeros(mocap.frames, 3);
for i = 1 : size(endeff.LEDxyz, 3)
   total = total + endeff.LEDxyz(:,:,i);
   pos = total / i;
end

% orientation of the end effector
orient = zeros(mocap.frames, 3);
for i = 1 : mocap.frames
    % find z-axis (just at one point for now for practice)
    A = endeff.LEDxyz(i,:,1);
    for j = 2 : size(endeff.LEDxyz, 3)
        A = [A ; endeff.LEDxyz(i,:,j)];    % build A one row at a time
    end
    checkA = isnan(A);
    if any(checkA(:) == 1)
        orient(i,:) = [NaN, NaN, NaN];      % if no sensor data, return NaN
    else
        [U,S,V] = svd(A);
        z = V'\[0 0 1]';
        led = endeff.LEDxyz(i,:,1);    % pick one of the LEDs, i.e. the "first" one
        xdir = (led - pos(i,:))';
        y = cross(z,xdir) / norm( cross(z,xdir) );
        x = cross(y,z) / norm( cross(y,z) );
        Rb2e = [x, y, z];
        Re2b = Rb2e';
        orient(i,:) = rotm2eul(Re2b);
    end
end

% state vector at each sample
endeff.x = [pos, orient];
endeff.t = mocap.t;


end





% %% CHECK CODE by plotting location of all LEDs at initial point in end effector coordinates
% for j = 1:3
%     poop = Hm2e*[mocap.topxyz(1,:,j),  1]';
%     pee(:,j) = poop(1:3);
% end
% for j = 1:5
%     fart = Hm2e*[mocap.effxyz(1,:,j),  1]';
%     butt(:,j) = fart(1:3);
% end
% 
% figure
% hold on
% plot3(pee(1,:), pee(2,:), pee(3,:), '*')
% plot3(butt(1,:), butt(2,:), butt(3,:), '*')
% hold off