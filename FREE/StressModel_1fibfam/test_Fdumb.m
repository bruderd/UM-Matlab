function F = test_Fdumb(r, params)

[x,y,z] = deal(r(1), r(2), r(3));
[c1, c2, c3] = deal(params.c(1), params.c(2), params.c(3));

F = c1*y*x^3 + c2*z*y^2 + c3*x*z;

end