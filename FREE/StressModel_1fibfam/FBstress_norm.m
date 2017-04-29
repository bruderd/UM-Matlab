% FBstress_norm
%   Finds the norm of the force balance equations at a point, so that
%   simulated annealing can use it as an objective function.

function Fnorm = FBstress_norm(x, u, params)

[F,J] = FBstress(x,u,params);
Fnorm = norm(F);

end