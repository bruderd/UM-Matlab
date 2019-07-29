function [ x , x_cm ] = alpha2x(alpha, params)
%alpha2x: Converts relative joint angles (alpha) to coordinates of joints (x)
%   and the coordinates of the links' centers of mass (x_cm).
%   x = [ x_0 ; y_0 ; x_1 ; y_1 ; ... ]

x = zeros( ( params.Nlinks + 1 ) * 2 ,  1 );
x_cm = zeros( params.Nlinks * 2 , 1 );

% if input is symbolic, so should output
if isa( alpha , 'sym' )
    x = sym(x);
    x_cm = sym(x_cm);
end

% convert to absolute joint angles (wrt vertical)
theta = alpha2theta(alpha);

% convert to coordinates of each joint (note there is 1 more joint than link)
for i = 1 : length(alpha)
   xim1 = x(2*(i-1)+1 : 2*i, 1); 
   x_cm(2*(i-1)+1 : 2*i, 1) = xim1 + params.l/2 * [ sin( theta(i) ) ; cos( theta(i) ) ];
   x(2*i+1 : 2*(i+1), 1) = xim1 + params.l * [ sin( theta(i) ) ; cos( theta(i) ) ];
end

end

