% troubleshooting with Ram

syms w

s = 1.0e+02 * -0.000208691882270;

[L_l, B_l, N_l] = deal(params.L_l, params.B_l, params.N_l); 
[L_r, B_r, N_r] = deal(params.L_r, params.B_r, params.N_r);

J_l = [pi*(B_l^2 - 3*(L_l+s)^2) / (2*pi*N_l+w)^2,...
    2*pi*(L_l+s)*((L_l+s)^2 - B_l^2) / (2*pi*N_l+w)^3];
    
J_r = [-pi*(B_r^2 - 3*(L_r-s)^2) / (2*pi*N_r-w)^2,...
   -2*pi*(L_r-s)*((L_r-s)^2 - B_r^2) / (2*pi*N_r-w)^3];

% call ezplot here to look at stuff...


%% Check the balance of forces of McKibbon muscle

[L, B, N] = deal(L_l, B_l, N_l); 
P1 = 5000;
P2 = 1000;

% use w in place of s since it's already a symbolic variable
poop = (B^2 - 3*(L-w)^2) * P1 - (B^2 - 3*(L+w)^2) * P2;

ezplot(poop)

