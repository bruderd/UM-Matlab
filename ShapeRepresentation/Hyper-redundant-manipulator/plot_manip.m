 function plot_manip( alpha, params )
%plot_manip: Displays what manipulator looks like given alpha
%   Detailed explanation goes here

[ X , Xcm ] = alpha2x(alpha, params);

x = [0; X(1:2:end)];
y = [0; X(2:2:end)];

figure
plot(x, y, '-o')
axis([-params.L, params.L, -0.5*params.L, 1.5*params.L])
set(gca,'Ydir','reverse')
xlabel('x(m)')
ylabel('y(m)')

end