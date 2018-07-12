% example: Spring-mass-damper system
clear

%% set parameters
params = struct;

params.n = 2;
params.maxDegree = 3;

params.m = 1;
params.b = 1;
params.k = 1;

% p is symbolic polynomial, psi is matrix of exponents: columns=monomial term, rows=dimension of x
[params.p, params.psi, params.N, xsym] = def_polyLift(params.n, params.maxDegree);

params.Ts = 0.01;
params.numTrials = 20;
params.x0max = [1, 1]';
params.tf = 10;

%% Simulate and find Koopman operator from "measurements"

[x,y] = get_snapshotPairs_SMD(params);

U = get_Koopman(x,y);

%% Calculate the infiniesimal generator as funtion of coeffients, and from data
Ldata = get_Ldata(U, params);
vecLdata = Ldata(:);

vecstackL = zeros(params.N^2, params.N*params.n);
for k = 1:params.N
    for j = 1:params.n
        Lkj = get_Lkj(k,j,params);
        
        % convert all the Lkj's into vectors and stack them horizontally
        vecLkj = Lkj(:);
        vecstackL(:, (k-1)*params.n + j) = vecLkj;
    end
end

%% solve for the coefficients, i.e. Eq. (18) from Mauroy and Gonclaves

W = pinv(vecstackL) * vecLdata;

% matrix of coefficents of monomials
w = reshape(W, [params.n, params.N]);

% dynamics (gives symbolic expression in terms of state)
vf = w * params.p; 
matlabFunction(vf, 'File', 'dynamics_sysid', 'Vars', {xsym});

%% Run simulatio of sysId'd system and compare results to real system

tspan = [0, params.tf];
x0sim = params.x0max*rand;
[tsysid, xsysid] = ode45(@(t,x) dynamics_sysid(x), tspan, x0sim);
[treal, xreal] = sim_SpringMassDamper(x0sim, params.tf, params.m, params.b, params.k);

% plot the results
figure
subplot(2,1,1)
plot(tsysid,xsysid)
title('Identified system')
subplot(2,1,2)
plot(treal,xreal)
title('Real system')

