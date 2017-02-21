function [F] = tensioneq_2fib_v2(x,u,x0)

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

    % Case 1: gama, betta > 0
    inputeq = u - P;
    force_balance = P*pi*r^2 - 2*(T_gama*cos(gama) + T_betta*cos(betta)) + F_elast;
    torque_balance = 2*r*(T_gama*sin(gama) + T_betta*sin(betta)) + M_elast;             % put (r) in front of tensions to fix units (2/2/2017)
    geometry_constraint1 = L/cos(gama) + r*(theta_gama0 + phi)/sin(gama);
    geometry_constraint2 = L/cos(betta) + r*(theta_betta0 + phi)/sin(betta);
    geometry_constraint3 = (theta_gama - theta_gama0) - phi;
    geometry_constraint4 = (theta_betta - theta_betta0) - phi;
    extra_constraint1 = P*r - (T_gama*sin(abs(gama)) + T_betta*sin(abs(betta)));
    

    F = [inputeq; force_balance; torque_balance; geometry_constraint1; geometry_constraint2; geometry_constraint3; geometry_constraint4; extra_constraint1];