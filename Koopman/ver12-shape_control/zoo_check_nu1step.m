% zoo_check_nu1step

function nu = zoo_check_nu1step( zeta0 , ydes , sys )

z0 = zeta0; % if zeta0 is already lifted use this and comment out line below
% z0 = sys.lift.full( zeta0 );

% % just for when ydes is actually zetades
% ydes = sys.lift.full( ydes );   

% % just for when ydes is state only but zeta includes 1 delay
% zetades = [ ydes ; zeta0( 2 * sys.params.n + 1 : 2 * sys.params.n + sys.params.m ) ];
% ydes = sys.lift.full( zetades );   

% C = sys.model.C * sys.model.B;
% d = -sys.model.C * sys.model.A * z0 + ydes;
C = eye( sys.params.mnu );
d = -sys.model.Beta * sys.model.A * z0 + sys.model.Beta * ydes;

nu = lsqlin( C , d , [] , [] );

end