function VF = EulerLagrange(L,X,Q_i,Q_e,R,par,varargin)
% EulerLagrange  Derives the system differential equations based on the
%                system Lagrangian L = T - V, where T is the kinetic and V
%                is the potential energy of the system and creates the
%                corresponding MATLAB function and Simulink block.
%
%   VF = EulerLagrange(L,X,Q_i,Q_e,R,par,varargin) uses the Lagrangian of a
%   system to derive its corresponding differential equations. The  
%   Lagrangian L needs to be defined symbolically in terms of the generalized 
%   coordinates and velocities X and system parameters par. Coordinates and
%   velocities need to be arranged in the form X = {q1 dq1 q2 dq2 ...}.
%   Q are the vectors of generalized forces. The Q_i are internal forces
%   and are functions of the generalized coordinates. Use the Q_e to 
%   describe external generalized forces that do not depend on the 
%   generalized coordinates but rather a yet to be determined control force
%   or torque, for example. R is the dissipation function as suggested by
%   Rayleigh. 
%
%   The function returns the vector field description of the derived
%   differential equations, VF, and the associated MATLAB function and 
%   Simulink block. The output is controlled by varargin, a cell array of 
%   strings:
%      'm':      create MATLAB function
%      's':      create Simulink block
%      filename: file name(s) for created output
%
%   Compatibility: This function has been successfully tested with the 
%   following releases
%      - R2013a
%      - R2014a through R2016b
%
% Example 1 - Projectile motion:
% syms g m x dx y dy z dz
% L   = m*(dx^2 + dy^2 + dz^2)/2 - m*g*z;
% X   = {x dx y dy z dz};
% Q_i = {0 0 0}; Q_e = {0 0 0};
% R   = 0;
% par = {g m};
% VF  = EulerLagrange(L,X,Q_i,Q_e,R,par);
%
% Example 2 - Rotary inverted pendulum:
% syms I m l Rl g al dal be dbe tau
% T   = I*dal^2/2 + m/2*(Rl^2*dal^2 + l^2*dal^2*sin(be)^2/3 + ...
%       l^2*dbe^2/4 - Rl*l*dal*dbe*cos(be)) + m*l^2*dbe^2/24;
% V   = m*g*l/2*cos(be);
% L   = T - V;
% X   = {al dal be dbe};
% Q_i = {-0.001*dal+0.001 -0.0001*dbe};
% Q_e = {tau  0};
% R   = 0;
% par = {I m l Rl g};
% VF  = EulerLagrange(L,X,Q_i,Q_e,R,par,'m','s','RotaryPendulum_ODE');
%
% Example 3 - Pendulum:
% syms th thd g m l k
% L   = m*l^2*thd^2/2 + m*g*l*cos(th);
% X   = {th thd};
% Q_i = {0}; Q_e = {0};
% R   = k*thd^2/2;
% par = {g m l k};
% VF  = EulerLagrange(L,X,Q_i,Q_e,R,par);
%
% Example 4 - Double pendulum:
% syms th1 th1d th2 th2d g m1 l1 m2 l2
% L   = (m1*l1^2*th1d^2 + ...
%        m2*(l1^2*th1d^2 + l2^2*th2d^2 + 2*l1*l2*th1d*th2d*cos(th1 - th2)))/2 + ...
%       (m1 + m2)*g*l1*cos(th1) + m2*g*l2*cos(th2);
% X   = {th1 th1d th2 th2d};
% Q_i = {0 0}; Q_e = {0 0};
% R   = 0;
% par = {g m1 l1 m2 l2};
% VF  = EulerLagrange(L,X,Q_i,Q_e,R,par,'m','s');
%
% Example 5 - RLC circuit:
% syms q qd Ri L C
% L   = L*qd^2/2 - q^2/(2*C);
% X   = {q qd};
% Q_i = {0}; Q_e = {0};
% R   = Ri*qd^2/2;
% par = {Ri L C};
% VF  = EulerLagrange(L,X,Q_i,Q_e,R,par);

% Copyright 2015-2016 The MathWorks, Inc.

% Define time as symbolic variable
syms t
rel_16a = '2016a';
% Reshape coordinates/velocities state vector to help with substitution
X2n     = reshape(X,2,[]);
numcoor = size(X2n,2);
% Create corresponding cell array of coordinates/velocities symfuns
rels = sort({version('-release'),rel_16a});
if strcmp(rels{1}, rel_16a)
    qt = arrayfun(@(ii) symfun(['q' int2str(ii) '(t)'],t),1:numcoor,'UniformOutput',false);
else
    qt = arrayfun(@(ii) symfun(sym(['q' int2str(ii) '(t)']),t),1:numcoor,'UniformOutput',false);
end
Xt2n = [qt; arrayfun(@(ii) diff(qt{1,ii},t),1:numcoor,'UniformOutput',false)];

R  = sym(R);  % take care of numeric zero
DE = [];
% Apply the Euler-Lagrange equation looping through all generalized
% coordinates and velocities
for ii = 1:numcoor
    % Differentiate w.r.t. coordinate/velocity
    L_q        = diff(L,X2n{1,ii});
    L_qdot     = diff(L,X2n{2,ii});
    R_qdot     = diff(R,X2n{2,ii});
    % Substitute coordinates/velocities by corresponding symfuns
    L_qt(t)    = subs(L_q,   X2n,Xt2n); %#ok<*AGROW>
    L_qtdot(t) = subs(L_qdot,X2n,Xt2n);
    R_qtdot(t) = subs(R_qdot,X2n,Xt2n);
    % Differentiate "velocity" term w.r.t. time
    L_qdotdt   = diff(L_qtdot,t);
    % Collect all Lagrange terms
    DE         = [DE; simplify(L_qdotdt - L_qt + R_qtdot) == ...
                                  subs(Q_i{ii} + Q_e{ii},X2n,Xt2n)];
end % of for

% Create vector field
[VF,Ysubs] = odeToVectorField(DE);

% Note: odeToVectorField does not preserve the order of the generalized
% coordinates (see documentation). Therefore,...
% Sort state vector to determine how re-ordering was done
[~,ind] = sort(Ysubs);
% Find indices of generalized coordinates and velocities
indq    = ind(numcoor+1:2*numcoor);
indqqd  = [indq'; indq'+1];
% Un-do re-ordering of generalized coordinates and velocities
exp_str = ['[Y[kk] $ kk = 1..',num2str(2*numcoor),']'];
Y       = evalin(symengine,exp_str);
VF      = subs(VF,Y,Y(indqqd(:)));
% Un-do re-ordering of vector field equations
VF      = VF(indqqd(:));

% Create corresponding MATLAB function and Simulink block
if (~isempty(varargin))
    if (~isempty(setdiff(varargin,{'m' 's'})))
        sysName = char(setdiff(varargin,{'m' 's'}));
    else
        sysName = 'my_ODE';
    end
    m_fileName  = strcat(sysName,'.m');
    s_fileName  = strcat(sysName,'_block');
    s_blockName = strcat(s_fileName,'/',sysName);
    for jj = 1:numel(varargin)
        Q_vars = Q_e(cellfun('isclass', Q_e, 'sym'));
        Q_vars = subs(Q_vars,X2n,Xt2n);
        switch varargin{jj}
            case {'m'} % create MATLAB function
                matlabFunction(VF,...
                    'file',     m_fileName,...
                    'vars',   [{'t','Y'}, Q_vars, par],...
                    'outputs', {'dYdt'});
            case {'s'} % create Simulink block
                new_system  (s_fileName)
                matlabFunctionBlock(s_blockName,VF,...
                    'functionName', sysName,...
                    'vars',       [{'t','Y'}, Q_vars, par],...
                    'outputs',     {'dYdt'});
                save_system (s_fileName)
                close_system(s_fileName)
        end % of switch
    end % of for
end % of if

% Express vector field in terms of original var names
VF = subs(VF,Y,X);
end % of function EulerLagrange