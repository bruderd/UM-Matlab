function product = complex_mult(z1,z2)
%complex_mult: Multiply two complex numbers specified as vectors
%   Detailed explanation goes here

real = z1(1) * z1(1) - z1(2) * z2(1);
im = z1(1) * z2(2) + z1(2) * z2(1);

product = [ real , im ];

end

