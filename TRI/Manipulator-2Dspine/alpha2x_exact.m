function x = alpha2x_exact( alpha, params )
%Converts alpha to x
%

p = params.p;
dl = params.dl;

x = zeros(3*p,1);

x(1:3, 1) = [0, dl, 0];
for i = 2:p
   xim1 = x(3*(i-2)+1 : 3*(i-1), 1); 
   theta = xim1(3,1) + alpha(i);
   x(3*(i-1)+1 : 3*i, 1) = xim1 + [sin(theta)*dl; cos(theta)*dl; alpha(i)];
end

end