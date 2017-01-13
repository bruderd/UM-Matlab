%% Explicit Dyanamics

function [f] = explicit_dynamics(x, u, params)

    n = params.n;
    m = params.m;
    
    [c1, c2, c3] = deal(params.Felast_consts(1), params.Felast_consts(2), params.Felast_consts(3));
    [c4, c5, c6] = deal(params.Melast_consts(1), params.Melast_consts(2), params.Melast_consts(3));
    
    [P0, gama0, r0, L0] = deal(params.x_rest(1), params.x_rest(2), params.x_rest(3), params.x_rest(4));
    [P, gama, r, L] = deal(x(1), x(2), x(3), x(4));
    

    
    dP = 0.5*(u - P);  % The constant in front of (u-P) is arbitrary
    dgama = -(2*P*dP*r^6*pi^2*cos(gama0)*cot(gama) - 4*P*dP*r^6*pi^2*cos(gama0)*cot(gama)^3 + pi*L*dP*r^2*cos(gama0)*tan(gama) - 2*pi*L*dP*r^2*cos(gama0)*cot(gama)^2*tan(gama))/(4*P^2*r^6*pi^2*cos(gama0) - L*L0*sin(gama)*tan(gama) + 20*P^2*r^6*pi^2*cos(gama0)*cot(gama)^2 + 16*P^2*r^6*pi^2*cos(gama0)*cot(gama)^4 + 2*pi*L*P*r^2*cos(gama0) - 4*pi*L*P*r^2*cos(gama0)*cot(gama)^2 + 2*pi*L*P*r^2*cos(gama0)*tan(gama)^2 - 6*pi*L0*P*r^4*cot(gama)*sin(gama) - 2*pi*L0*P*r^2*sin(gama)*tan(gama) - 4*pi*L*P*r^2*cos(gama0)*cot(gama)^2*tan(gama)^2 + 4*pi*L*P*r^2*cos(gama0)*cot(gama)*tan(gama) + 4*pi*L*P*r^2*cos(gama0)*cot(gama)^3*tan(gama) + 4*pi*L0*P*r^2*cot(gama)^2*sin(gama)*tan(gama));
    dr = -(r^2*(2*P*dP*r^5*pi^2*cos(gama0) + pi*L*dP*r*cos(gama0) - pi*L0*dP*r*sin(gama)*tan(gama) + 6*P*dP*r^5*pi^2*cos(gama0)*cot(gama)^2 + 4*P*dP*r^5*pi^2*cos(gama0)*cot(gama)^4 - 2*pi*L*dP*r*cos(gama0)*cot(gama)^2 - 2*pi*L0*dP*r^3*cot(gama)*sin(gama) + pi*L*dP*r*cos(gama0)*tan(gama)^2 - 2*pi*L*dP*r*cos(gama0)*cot(gama)^2*tan(gama)^2 + 2*pi*L0*dP*r*cot(gama)^2*sin(gama)*tan(gama)))/(4*P^2*r^6*pi^2*cos(gama0) - L*L0*sin(gama)*tan(gama) + 20*P^2*r^6*pi^2*cos(gama0)*cot(gama)^2 + 16*P^2*r^6*pi^2*cos(gama0)*cot(gama)^4 + 2*pi*L*P*r^2*cos(gama0) - 4*pi*L*P*r^2*cos(gama0)*cot(gama)^2 + 2*pi*L*P*r^2*cos(gama0)*tan(gama)^2 - 6*pi*L0*P*r^4*cot(gama)*sin(gama) - 2*pi*L0*P*r^2*sin(gama)*tan(gama) - 4*pi*L*P*r^2*cos(gama0)*cot(gama)^2*tan(gama)^2 + 4*pi*L*P*r^2*cos(gama0)*cot(gama)*tan(gama) + 4*pi*L*P*r^2*cos(gama0)*cot(gama)^3*tan(gama) + 4*pi*L0*P*r^2*cot(gama)^2*sin(gama)*tan(gama));
    dL = (L0*r^2*pi*sin(gama)*(L*dP*tan(gama) - 2*L*dP*cot(gama)^2*tan(gama) + 2*pi*P*dP*r^4*cot(gama) - 4*pi*P*dP*r^4*cot(gama)^3))/(4*P^2*r^6*pi^2*cos(gama0) - L*L0*sin(gama)*tan(gama) + 20*P^2*r^6*pi^2*cos(gama0)*cot(gama)^2 + 16*P^2*r^6*pi^2*cos(gama0)*cot(gama)^4 + 2*pi*L*P*r^2*cos(gama0) - 4*pi*L*P*r^2*cos(gama0)*cot(gama)^2 + 2*pi*L*P*r^2*cos(gama0)*tan(gama)^2 - 6*pi*L0*P*r^4*cot(gama)*sin(gama) - 2*pi*L0*P*r^2*sin(gama)*tan(gama) - 4*pi*L*P*r^2*cos(gama0)*cot(gama)^2*tan(gama)^2 + 4*pi*L*P*r^2*cos(gama0)*cot(gama)*tan(gama) + 4*pi*L*P*r^2*cos(gama0)*cot(gama)^3*tan(gama) + 4*pi*L0*P*r^2*cot(gama)^2*sin(gama)*tan(gama));

    f = [dP; dgama; dr; dL];     


         