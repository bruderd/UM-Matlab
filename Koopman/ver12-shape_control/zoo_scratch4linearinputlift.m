% zoo_scratch4linearinputlift

% Put something like this in the basis definition function
% NEW: multiply input by monomials of state, but no powers of input
if obj.liftinput == 1
    basis_diag = kron( eye(obj.params.m) , fullBasis );
    fullBasis = [ fullBasis ; basis_diag * u ];
    zeta = [ zeta ; u ]; % I think this is needed to make it run
end