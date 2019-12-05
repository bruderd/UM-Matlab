% zoo_check_nu1step

function nu = zoo_check_nu1step_2( zx , ydes , sys )

H = sys.model.C * sys.model.Beta_pinv;
f = ydes - sys.model.C * sys.model.A * zx - sys.model.C * sys.traindata.emean';

nu = lsqlin( H , f , [] , [] );

end