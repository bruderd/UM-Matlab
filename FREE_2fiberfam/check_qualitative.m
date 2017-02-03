% Check qualitative behavior of FREEs

Num = 5;
T = 1;
qual = [0 0 0 0];

for i = -17:17
    for j = -17:17
        L0 = 5;
        gama0 = 5*i + 1;
        betta0 = 5*j + 2;
        x_rest = [0.0001, deg2rad(gama0), deg2rad(betta0), 3/16, L0, 0]';
        [x_star, u_star] = main_2fiberfam_func(x_rest, Num, T);
        
        close all
        
        qual = [qual; gama0, betta0, sign(x_star(5,end)-L0), sign(x_star(6,end))];
    end
end