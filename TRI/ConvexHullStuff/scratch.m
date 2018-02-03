clear;

M = [eye(4);...
    1 1 0 0; 0 1 1 0; 0 0 1 1; 1 0 0 1; 1 0 1 0; 0 1 0 1;...
    0 1 1 1; 1 0 1 1; 1 1 0 1; 1 1 1 0];

x1 = [-5,2];
x2 = [-4,-4];
x3 = [2,6];
x4 = [2,-4];

X = [x1; x2; x3; x4];

V = M*X;

vx = V(:,1);
vy = V(:,2);

k = convhull(vx,vy);

figure
hold on
plot(X(:,1),X(:,2),'r*');
plot(vx(k), vy(k), 'b-');
hold off