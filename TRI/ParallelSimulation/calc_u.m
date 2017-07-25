function u = calc_u( t )
%Calculate the input to the system, which is a vector of pressures
%   Detailed explanation goes here

% u = [1 0 1 0]' * 1000*(sin(t) + 1);

% u = [1 0 1 0]' * 1000 * (sin(pi*t/5) + 1);

u = [0.6 0 0.4 0.8]' * 1000 * (sin(pi*t/5)*exp(-t*1e-2) + 1);

end

