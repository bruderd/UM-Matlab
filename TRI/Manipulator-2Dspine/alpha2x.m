function x = alpha2x( alpha, params )
%Converts alpha to x
%   Detailed explanation goes here

p = params.p;
dl = params.dl;

x = zeros(3*p,1);

x(1:3, 1) = [0, dl, 0];
for i = 2:p
   xim1 = x(3*(i-2)+1 : 3*(i-1), 1); 
   theta = xim1(3,1) + alpha(i);
   x(3*(i-1)+1 : 3*i, 1) = xim1 + [-theta*dl; (1-0.5*theta^2)*dl; alpha(i)];
end

end

