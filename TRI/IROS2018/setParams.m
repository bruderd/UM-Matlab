function params = setParams(Gama, R, L, d, p, Pmax)
%setParams: Creates a struct containing all parameters of the problem
%   Detailed explanation goes here

% if size(Gama,2) ~= size(R,2) ~= size(L,2) ~= size(d,2) ~= size(p,2)
%     disp('Column number of all inputs must equal the number of FREEs in parallel')
%     return
% end

params = struct;

params.Gama = Gama;
params.R = R;
params.L = L;
params.d = d;
params.p = p;
params.Pmax = Pmax;

params.B = abs(params.L ./ cos(params.Gama));   % fiber length (must be positive))
params.N = -params.L ./ (2*pi*params.R) .* tan(params.Gama); % total fiber windings in revolutions (when relaxed)

params.num = size(Gama,2);  % total number of FREEs

end

