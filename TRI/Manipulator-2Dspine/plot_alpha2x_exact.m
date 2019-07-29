 function plot_alpha2x_exact( alpha, params )
%plot_alpha2x: Displays what manipulator looks like given alpha
%   Detailed explanation goes here

X = alpha2x_exact(alpha, params);

x = [0; X(1:3:end)];
y = [0; X(2:3:end)];
theta = X(3:3:end);

figure
plot(x, y, '-o')
axis([-params.L, params.L, 0, params.L])
set(gca,'Ydir','reverse')
xlabel('x(m)')
ylabel('y(m)')

end