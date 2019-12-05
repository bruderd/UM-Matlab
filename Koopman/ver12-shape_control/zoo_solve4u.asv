function u = zoo_solve4u( nu , x , u0 , sys )
%zoo_solve4u: Solves for input u given pseudoinput nu
%   Detailed explanation goes here

e = sys.model.Beta_pinv * nu + sys.traindata.emean';
Fxu = e - sys.model.B * sys.basis.zu;
Fu = subs( Fxu , sys.params.x , x );

% sum all of the equations (note, only works for 1-dim u) 
Fu_sum = sum( Fu );

% convert to polynomial so we can solve using the roots function
Fu_poly = fliplr( coeffs( Fu_sum , sys.params.u ) );

% solve for the input, make sure it's positive
Fcost = matlabFunction( Fu' * Fu );  % want Fu = 0 
u = fmincon( Fcost , u0 , [] , [] , [] , [] , -Inf , Inf );

% % solve for u
% u = roots( Fu_poly );

end

