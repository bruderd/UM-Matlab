function params = setInvKin(params)
%setInvKin: Creates function to evaluate inverse kinematics (x2q) based
%on user supplied definition.
%   Detailed explanation goes here

num = params.num;
d = params.d;
L = params.L;

%% (USER EDIT) Symbolically defined inverse kinematics function

q = sym('q', [2*params.num, 1], 'real'); 
x = sym('x', [6,1], 'real');

% % 2 DOF rotation and translation rig (axes aligned)
% q(1:2:end) = (x(3) - params.deff(3))*ones(params.num,1);  % dl
% q(2:2:end) = x(4)*ones(params.num,1);  % dphi
% 
% % 2 DOF rotation and translation rig (axes not aligned)
% for i = 1 : num
%     q(2*(i-1) + 1) = -L(i) + sqrt( (L(i)+x(3))^2 + ( x(4)*(d(1,i) + d(2,i)) )^2);
%     q(2*(i-1) + 2) = x(4);
% end

% 6 DOF module with compliant spine
for i = 1 : num
    q(2*(i-1) + 1) = -L(i) + sqrt( ( L(i) - sin(x(5))*d(1,i) + cos(x(5))*sin(x(6))*d(2,i) )^2 + (x(4)*(d(1,i) + d(2,i)))^2);
    q(2*(i-1) + 2) = x(4);
end

%% (USER EDIT) Symbolically definded coupling between states

% define which components are independent (all others can be written as a function of these components)
params.xindselect = [0 0 0 1 1 1];  % specify which states are independent
params.xindnum = length( params.xindselect(params.xindselect == 1) );   % number of independent states

% % 2 DOF rotation and translation rig (axes aligned)
% xind = sym('xind', [2,1], 'real'); % size should be number of independent components
% xcoupled = [0,...
%             0,...
%             xind(1),...
%             xind(2),...
%             0,...
%             0,];


% % 2 DOF rotation and translation rig (axes not aligned)
% xind = sym('xind', [2,1], 'real'); % size should be number of independent components
% xcoupled = [0,...
%             0,...
%             -L(1) + sqrt( (L(1) + xind(1))^2 - (xind(2)*(d(1,i) + d(2,i)))^2 ),...
%             xind(2),...
%             0,...
%             0,];
            
        
% 6 DOF module with compliant spine
xind = sym('xind', [3,1], 'real'); % size should be number of independent components
xcoupled = [L(1)/2 * (cos(xind(1))*sin(xind(2))*cos(xind(3)) + sin(xind(1))*sin(xind(3))),...
            L(1)/2 * (sin(xind(1))*sin(xind(2))*cos(xind(3)) - cos(xind(1))*sin(xind(3))),...
            L(1)/2 * (cos(xind(2))*cos(xind(3)) - 1),...
            xind(1),...
            xind(2),...
            xind(3)];

%% Convert symbolic expression to matlab function
matlabFunction(q, 'File', 'x2q', 'Vars', {x});  
matlabFunction(xcoupled, 'File', 'xind2xcoupled', 'Vars', {xind});

end

