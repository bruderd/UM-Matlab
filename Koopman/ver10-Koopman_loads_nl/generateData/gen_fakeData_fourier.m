% gen_fakeData_fourier

tspan = 0:0.01:100;
x0 = zeros(3,1);
[tsol, xsol] = ode45(@(t,x) vf(x), tspan, x0);

function xdot = vf(x)

xdot = [sin(x(2)); cos(x(1)); -sin(x(2))]; 

end