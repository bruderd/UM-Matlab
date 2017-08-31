function u = calc_u( t )
%Calculate the input to the system, which is a vector of pressures
%   Detailed explanation goes here

% u = [1 0 0 1]' * 1e6;

% u = [1 0 1 0]' * 1000*(sin(t) + 1);

% u = [1 0 1 0]' * 1000 * (sin(pi*t/5) + 1);

% u = [0.6 0 0.4 0.8]' * 1e8 * (sin(pi*t/5)*exp(-t*2e-1) + 1);

u = [0.6 0 0.4 0.8]' * 1e8 * (1 - exp(-t * 2e0));
end

