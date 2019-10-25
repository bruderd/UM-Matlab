% zoo_check_e1step

function e = zoo_check_e1step( zeta0 , ydes , sys )

z0 = zeta0; % if zeta0 is already lifted use this and comment out line below
% z0 = sys.lift.full( zeta0 );
% ydes = sys.lift.full( ydes );   % just for when ydes is actually zetades

% C = sys.model.C;
% d = -sys.model.C * sys.model.A * z0 + ydes;

% for when everything is alrady lifted to the zx space
C = eye( sys.params.nzx );
d = -sys.model.A * z0 + ydes;

e = lsqlin( C , d , [] , [] );

end