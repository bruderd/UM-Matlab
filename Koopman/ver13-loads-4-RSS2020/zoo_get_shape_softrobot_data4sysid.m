% zoo_get_shape_softrobot_data4sysid
%
% Transforms the mocap data from the soft robot test rig into shape
% polynomial coefficients.
% Use data file: 

% points along the body
s = [1:3]' / 3;

% A matrix (assuming 3rd order polynomials for shape)
A = [ s , s.^2 , s.^3 ];
Adagger = pinv(A);

for i = 1 : length(train)
    
    % isolate xyz components of the markers
    t.x = train{i}.y(:,1:3:end);
    t.y = train{i}.y(:,2:3:end);
    t.z = train{i}.y(:,3:3:end);
    v.x = val{i}.y(:,1:3:end);
    v.y = val{i}.y(:,2:3:end);
    v.z = val{i}.y(:,3:3:end);
    
    % define the polynomial coefficients at each timestep for each coordinate's polynomial
    px_train = t.x * Adagger';
    py_train = t.y * Adagger';
    pz_train = t.z * Adagger';
    px_val = v.x * Adagger';
    py_val = v.y * Adagger';
    pz_val = v.z * Adagger';
    
    % stack these back together into a single data array
    coeffs_train = [ px_train , py_train , pz_train ];
    coeffs_val = [ px_val , py_val , pz_val ];
    
    % create a training data struct and validation data struct
    train{i}.y = coeffs_train;
    val{i}.y = coeffs_val;
    
end