% gen_fakeData_sin

tspan = 0:0.01:100;
u = sin(0.1 * tspan)';
x0 = 0;
[tsol, xsol] = ode45(@(t,x) vf(x, get_u(t, tspan, u)), tspan, x0);

function xdot = vf(x,u)

xdot = sin(2*pi*u); 

end

function unow = get_u(t, tspan, u)

unow = interp1(tspan, u, t)';

end