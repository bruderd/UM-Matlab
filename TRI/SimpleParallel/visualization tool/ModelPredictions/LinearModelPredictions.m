% LinearModelPredictions.m
%   Determines the displacements and corresponding forces for a 1-fiber FREE,
%   according to the volumetric Jacobian model (linear).
%
%   Must specify parameters in setParms.m before running this script.

% Give parameters local names
[Gama, R, L, B, N] = deal(params.GamaDeg, params.R, params.L, params.B, params.N);
[Kf, Km] = deal(params.kelast(1), params.kelast(2));
[F_load, M_load] = deal(params.load(1), params.load(2));
res = params.res;
Ptest = params.Ptest;

smin = -0.5*L;
smax = 0.5*L;
wmin = -2*pi;
wmax = 2*pi;
srange = linspace(smin,smax,res);
wrange = linspace(wmin,wmax,res);

[s,w] = meshgrid(srange, wrange);

% Calculate the volume
V = pi * (L+s) .* (B.^2 - (L+s).^2) ./ (2*pi*N + w).^2;

% Calculate the volume Jacobian
dVds = pi*(B^2 - 3.*(L+s).^2) ./ (2*pi*N + w).^2;
dVdw = 2*pi*((L+s) .* ((L+s).^2 - B^2)) ./ (2*pi*N + w).^3;

% Calculate elastomer forces
Felast = Kf * s;
Melast = Km * w;

% Multiply by pressure to get the forces
for i = 1:length(Ptest)
    Force(:,:,i) = Ptest(i)*10^3 * dVds + Felast;
    Moment(:,:,i) = Ptest(i)*10^3 * dVdw + Melast;
    
%     % No elastomer
%     Force(:,:,i) = Ptest(i)*10^3 * dVds;
%     Moment(:,:,i) = Ptest(i)*10^3 * dVdw;
end

%% Convert matrices to vectors
s_vec = matrix2vector(s);
w_vec = matrix2vector(w);
V_vec = matrix2vector(V);
F = matrix2vector(Force(:,:,1));
M = matrix2vector(Moment(:,:,1));
P = ones(length(s)^2,1) * Ptest(1);
if length(Ptest) > 1;
    for j = 2:length(Ptest)
        %     s_vec( (j-1)*length(s)+ 1 : j*length(s), 1) = matrix2vector(s);
        %
        %
        %     F( (j-1)*length(s)+ 1 : j*length(s), 1) = matrix2vector(Force(:,:,j));
        %     M( (j-1)*length(s)+ 1 : j*length(s), 1) = matrix2vector(Force(:,:,j));
        
        s_vec = [s_vec; matrix2vector(s)];
        w_vec = [w_vec; matrix2vector(w)];
        V_vec = [V_vec; matrix2vector(V)];
        F = [F; matrix2vector(Force(:,:,j))];
        M = [M; matrix2vector(Moment(:,:,j))];
        P = [P; ones(length(s)^2,1) * Ptest(j)];
    end
end

%% Creat csv from data
Matrix = [zeros(1, 3+length(Ptest));...
          s_vec, w_vec, F, M, P, V_vec, zeros(length(s_vec), length(Ptest) - 3)];

Matrix(1, 1 : 3+length(Ptest)) = [Gama, L, R, Ptest];     % add headers to top row

csvwrite('LinModelData.csv', Matrix);
