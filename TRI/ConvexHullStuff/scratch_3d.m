clear;

M = [eye(4);...
    1 1 0 0; 0 1 1 0; 0 0 1 1; 1 0 0 1; 1 0 1 0; 0 1 0 1;...
    0 1 1 1; 1 0 1 1; 1 1 0 1; 1 1 1 0];

x1 = [1,1,1];
x2 = [-1,-1,1];
x3 = [1,-1,-1];
x4 = [-1,1,-1];

X = [x1; x2; x3; x4];

V = M*X;

vx = V(:,1);
vy = V(:,2);
vz = V(:,3);

k = convhull(vx,vy,vz);

figure
hold on
plot3(X(:,1),X(:,2), X(:,3),'r*');
% plot3(vx(k), vy(k), vz(k), 'b-');
h = trisurf(k,vx,vy,vz,'Facecolor','red','FaceAlpha',0.1)
hold off