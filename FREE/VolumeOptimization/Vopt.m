function [ x_star, cost_star, c_star, ceq_star ] = Vopt( Gama, L, R )
%Vopt: Solves the FREE volume maximization problem.
%   Optimization results vary too much depending on initial condition to be
%   useful :(

params = setParams(Gama, L, R);

% x0 = [L, 0, R, 0, 0, 0, 0, 0, 0];
x0 = [L, 0, R, rand(1,6)];
lb = [L/4, -inf, R-(1e-6), -inf, -inf, -inf, -inf, -inf, -inf];
ub = [params.S, inf, inf, inf, inf, inf, inf, inf, inf];
options = optimoptions('fmincon', 'algorithm', 'interior-point', 'Display','iter');
options.MaxFunctionEvaluations = 6000;

x_star = fmincon(@(x) costwgrad(x),x0,[],[],[],[],lb,ub, @(x) constraints(x,params), options);
[c_star, ceq_star] = constraints( x_star, params);
cost_star = cost(x_star);

% Evauluate the radius polynomial
points = linspace(0,x_star(1),100);
poly = x_star(3:end);
r = fliplr(poly);
p_star =  polyval(r,points);

% Plot the radius polynomial
figure
plot(points, p_star);
% axis equal;
% xlim([0, 1.25*L]);
% ylim([0, 1.5*R])

figure
plot(points, p_star);
axis equal;
xlim([0, 1.25*L]);
ylim([0, 3*R])


end

