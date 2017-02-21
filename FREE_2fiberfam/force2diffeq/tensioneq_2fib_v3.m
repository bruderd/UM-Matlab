% PROBABLY BAD: does not produce the same tensions as v2
% v3: Uses force balance equations on arbitryry lengh of tube to calculate
%   tension forces (the old way).

function [F] = tensioneq_2fib_v3(x,u,x0)

    gama0 = x0(2);
    betta0 = x0(3);
    r0 = x0(4);
    L0 = x0(5);
    phi0 = x0(6);
    T_gama0 = x0(7);
    T_betta0 = x0(8);

    P = x(1);
    gama = x(2);
    betta = x(3);
    r = x(4);
    L = x(5);
    phi = x(6);
    T_gama = x(7);
    T_betta = x(8);
    
    %simple elastomer constants
    c1 = 7;
    c4 = 7;
    
    %Elastomer forces: simpler version
%     F_elast = c1*(L0-L);    
%     M_elast = c4 * phi;    % removed (-1) becasue I changed phi up there (1/26/2017)
    F_elast = 0;    
    M_elast = 0;    % removed (-1) becasue I changed phi up there (1/26/2017)
    

    theta_gama0 = -tan(gama0)*L0/r0;     % (-) fixes sign convention (1/28/2017)
    theta_betta0 = -tan(betta0)*L0/r0;   % (-) fixes sign convention (1/28/2017)
    theta_gama = -tan(gama)*L/r;       % (-) fixes sign convention (1/28/2017)
    theta_betta = -tan(betta)*L/r;       % (-) fixes sign convention (1/28/2017)

%     % force balance equations
    inputeq = u - P;
%     force_balance = P*pi*r^2 - 2*(T_gama*cos(gama) + T_betta*cos(betta)) + F_elast;
%     torque_balance = 2*r*(T_gama*sin(gama) + T_betta*sin(betta)) + M_elast;             % put (r) in front of tensions to fix units (2/2/2017)
%     geometry_constraint1 = L/cos(gama) + r*(theta_gama0 + phi)/sin(gama);
%     geometry_constraint2 = L/cos(betta) + r*(theta_betta0 + phi)/sin(betta);
%     geometry_constraint3 = (theta_gama - theta_gama0) - phi;
%     geometry_constraint4 = (theta_betta - theta_betta0) - phi;
%     extra_constraint1 = P*r - (T_gama*sin(abs(gama)) + T_betta*sin(abs(betta)));

    % tension equations
    nrat = floor(abs(tan(gama)/tan(betta)));
    psi = pi*(tan(abs(gama))/tan(abs(betta)) - nrat);
    tension1 = -4*pi*P*r^2*cot(abs(gama)) + 2*T_gama*sin(abs(gama)) + T_betta*sin(abs(betta))*(1 - cos(pi*tan(abs(betta))/tan(abs(gama))));
    tension2  = -4*pi*P*r^2*cot(abs(betta)) + 2*T_betta*sin(abs(betta)) + T_gama*sin(abs(gama))*((2*nrat+1) - cos(psi));
    

    F = [inputeq, gama - gama0, betta - betta0, r - r0, L - L0, phi - phi0, tension1, tension2];