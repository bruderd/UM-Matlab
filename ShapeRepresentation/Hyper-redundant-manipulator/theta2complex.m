function complex = theta2complex(theta)
%angle2complex: Converts an angle relative to z-axis to a point on the complex unit circle
%   Note that the answer is an array [a b] for the complex number a+ib

a = sin( theta );
b = cos( theta );

complex = [ a , b ];

end

