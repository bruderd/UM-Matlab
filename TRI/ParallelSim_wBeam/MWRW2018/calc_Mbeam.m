function tauBeam = calc_Mbeam( x_eul, params )
%calc_Mbeam: calculates the moment exerted on the end effector by central
%spine beam.
%   Detailed explanation goes here

[Ebeam, Ibeam] = deal(params.Ebeam, params.Ibeam);
x_cart = euler2cart(x_eul, params);
[x,y,z] = deal(x_cart(1), x_cart(2), x_cart(3) + params.Lspine);
R_10 = R10(x_eul);
R_01 = R01(x_eul);



d = sqrt(x^2 + y^2 + z^2);  % distance between eff and origin
p = [0 0 1]';       % normal z vector in 0-frame
q = R_10 * p;        % normal z vecotr in 1-frame
gama = acos(p' * q);    % angle between normal vectors in 0,1-frame
rho = (d/2) * ( 1/sin(gama/2) ); % radius of curvature of beam

Mbeam = Ebeam*Ibeam/rho;    % beam moment

% case to prevent beta = NaN
if x == 0 && y == 0
    beta = 1;   % value doesn't matter since it will be zeroed out
else
    beta = atan(y/x);
end

tauBeam0 = [-Mbeam * sin(beta); Mbeam * cos(beta); 0];  % beam moment in frame-0
tauBeam = R_01 * tauBeam0;     % beam moment in frame-1


end

