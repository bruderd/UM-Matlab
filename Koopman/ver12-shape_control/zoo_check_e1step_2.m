% zoo_check_nu1step

function e = zoo_check_e1step_2( zx , ydes , sys )

H = sys.model.C;
f = ydes - sys.model.C * sys.model.A * zx;

e = lsqlin( H , f , [] , [] );

end