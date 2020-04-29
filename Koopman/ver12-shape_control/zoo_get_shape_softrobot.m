% zoo_get_shape_softrobot
%
% Transforms the mocap data from the soft robot test rig into shape
% polynomial coefficients.
% Use data file: /v2_12hz/v3_3mods_rampNhold_0-10v_ramp1s_14

% isolate xyz components of the markers
dat.x = sysidData.Y(:,1:3:end);
dat.y = sysidData.Y(:,2:3:end);
dat.z = sysidData.Y(:,3:3:end);

% points along the body
s = [0:6]' / 6;

% A matrix (assuming 3rd order polynomials for shape)
A = [ s , s.^2 , s.^3 ];
Adagger = pinv(A);


% define the polynomial coefficients at each timestep for each coordinate's polynomial
px = dat.x * Adagger';
py = dat.y * Adagger';
pz = dat.z * Adagger';

% stack these back together into a single data array
coeffs = [ px , py , pz ];

% create a training data struct and validation data struct
% train.y = coeffs(150:600,:); 
train.y = sysidData.Y(150:600 , [7:9,13:15,19:21] );
train.u = sysidData.U(150:600,:);
train.t = sysidData.T(150:600,:) - sysidData.T(150,:);

% val.y = coeffs(601:end,:); 
val.y = sysidData.Y(601:end , [7:9,13:15,19:21]); 
val.u = sysidData.U(601:end,:);
val.t = sysidData.T(601:end,:) - sysidData.T(601,:);