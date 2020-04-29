% def_multilinear_lift
%
% Defines a multilinear lifting function.
% Valid only for SISO systems

function def_multilinear_lift( horizon )

% symbolic sequence of inputs
U = sym( 'U' , [ 1 , horizon ] );

exponents = fliplr( permn( [0,1] , horizon ) );

monomials = prod( U .^ exponents , 2);

% turn symbolic expression into matlab function (accepts row vectors as inputs)
matlabFunction( monomials , 'File' , 'multilinear_lift' , 'Vars' , {U} );

end