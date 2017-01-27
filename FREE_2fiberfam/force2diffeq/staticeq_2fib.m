function [F] = staticeq_2fib(x,u,x0)

    P = u;
    
    gama0 = x0(2);
    betta0 = x0(3);
    r0 = x0(4);
    L0 = x0(5);
    phi0 = x0(6);

    gama = x(2);
    betta = x(3);
    r = x(4);
    L = x(5);
    phi = x(6);
    
    %simple elastomer constants
    c1 = 7;
    c4 = 7;
    
    %simpler version
    F_elast = c1*(L0-L);    
    M_elast = c4 * phi;    % removed (-1) becasue I changed phi up there (1/26/2017)
    
    T_gama = (2*(pi*P*r^2*cot(gama) + P*r^2*pi*cos((pi*tan(betta))/tan(gama))*cot(betta)))/(sin(gama)*(cos((pi*tan(betta))/tan(gama)) + 1));
    T_betta = (2*(pi*P*r^2*cot(betta) - pi*P*r^2*cot(gama)))/(sin(betta) + cos((pi*tan(betta))/tan(gama))*sin(betta));
    theta_gama0 = tan(gama0)*L0/r0;
    theta_betta0 = tan(betta0)*L0/r0;

    
    inputeq = P - u;
    force_balance = P*pi*r^2 - (T_gama*cos(gama) + T_betta*cos(betta)) + F_elast;   
    torque_balance = (T_gama*sin(gama) + T_betta*sin(betta)) + M_elast;               
    geometry_constraint1 = L/cos(gama) - r*(theta_gama0 + phi)/sin(gama);
    geometry_constraint2 = L/cos(betta) - r*(theta_betta0 + phi)/sin(betta);
    geometry_constraint3 = (L/r)*tan(gama) - (L0/r0)*tan(gama0) - phi;
    geometry_constraint4 = (L/r)*tan(betta) - (L0/r0)*tan(betta0) - phi;
    
    F = [inputeq; force_balance; torque_balance; geometry_constraint1; geometry_constraint2; geometry_constraint3];