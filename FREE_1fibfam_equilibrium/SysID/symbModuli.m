% Symbolically solve for E, G, and T in terms of the states of the system
clear

syms P gama r L phi T G E P_rest gama_rest r_rest L_rest phi_rest t

%% With tube thickness included
HB = pi*P*r^2*cot(gama) - sin(gama)*T - pi*r*cot(gama)*((r-r_rest)/r_rest)*E*t;
FB = P*pi*r^2 - 2*cos(gama)*T - 2*pi*r*((L-L_rest)/L_rest)*E*t;
TB = sin(gama)*T - (pi*r^2*phi/L)*G*t;

solution = solve([HB; FB; TB], [T, E, G]);

%% Without tube thickness included
HB = pi*P*r^2*cot(gama) - sin(gama)*T - pi*r*cot(gama)*((r-r_rest)/r_rest)*E;
FB = P*pi*r^2 - 2*cos(gama)*T - 2*pi*r*((L-L_rest)/L_rest)*E;
TB = sin(gama)*T - (pi*r^2*phi/L)*G;

solution = solve([HB; FB; TB], [T, E, G]);

pole = solve((L*r_rest*sin(gama) - L_rest*r_rest*sin(gama) - L_rest*r*cos(gama)*cot(gama) + L_rest*r_rest*cos(gama)*cot(gama)), [r,gama,L]);