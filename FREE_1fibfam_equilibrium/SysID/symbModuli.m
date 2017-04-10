% Symbolically solve for E, G, and T in terms of the states of the system

syms P gama r L phi T G E P_rest gama_rest r_rest L_rest phi_rest t

HB = pi*P*r^2*cot(gama) - sin(gama)*T - pi*r*cot(gama)*((r-r_rest)/r_rest)*E*t;
FB = P*pi*r^2 - 2*cos(gama)*T - 2*pi*r*((L-L_rest)/L_rest)*E*t;
TB = sin(gama)*T - (pi*r^2*phi/L)*G*t;

solution = solve([HB; FB; TB], [T, E, G]);
