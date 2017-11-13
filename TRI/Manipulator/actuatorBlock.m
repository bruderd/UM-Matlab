function [Z, Vdot] = actuatorBlock( P, qdot, x, params)
%actuatorBlock: Relates inputs (pressure) to forces generated by the actuators (FREEs)
%   Currently using the linear FREE model. NOTE: This is only a
%   placeholder, this function can be used to implement any FREE model by
%   replacing the contents of this function.

p = params.p;   % total modules in manipulator
n = params.n;   % number of actuators in each module (a vector)
B = params.B;   % FREE fiber length (array)
N = params.Nf;   % FREE total fiber windings in revolutions (when relaxed) (array)
L = params.L;   % length of central spine of each module

% Define the extension and twist of each FREE with respect to x
q = x2q(x , params);
s = q(1:2:end);
w = q(2:2:end);

%% Define volume Jacobain
Jv = zeros(sum(n), 2*sum(n));      % initialize Jv
for i = 1:sum(n)
    Jv_ki = [(pi*(B(i)^2 - 3*(L+s(i))^2) / (2*pi*N(i)+w(i))^2),...
             2*pi*(L+s(i))*((L+s(i))^2 - B(i)^2) / (2*pi*N(i)+w(i))^3];
         
    Jv( i, 2*(i-1)+1 : 2*i ) = Jv_ki;     % stack volume jacobian for each actuator diagonally to form Jv.   
end

%% Define elastomer stiffness and damping matrices

Kc = zeros(length(q), length(q));   % elastomer stiffness matrix
Kd = zeros(length(q), length(q));   % elastomer damping matrix

% elastomer spring constants (could set these values in setParams, but I wanted to keep that script agnostic about the FREE model for now).
kelast = [-4e1, -4e-2]';     % [axial stiffness, rotational stiffness]'

% insert stiffnesses as diagonal elements of Kc (really janky way to do this...)
for j = 1:length(q)/2
    Kc(2*j-1:2*j, 2*j-1:2*j) = [kelast(1), 0; 0, kelast(2)];
end

%% Set output values
Vdot = Jv * qdot;
Z = Jv' * P + Kc * q + Kd * qdot;



end

