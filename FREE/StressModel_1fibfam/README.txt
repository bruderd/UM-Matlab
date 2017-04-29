README for StressModel_1fibfam
 
setParams.m - Sets the values of all of the FREE parameters as well as the measured values of P, dL, and dphi for several fiber angles.
Calls:
-importfile_loads.m
-FBsymbolic
    
importfile_loads.m - imports measured values of P, dL, and dphi from csv files stored in the ‘data’ subfolder.

FBsymbolic.m - Defines the set of force balance equations and its jacobian symbolically, then uses matlabFunction to turn the symbolic expressions into the functions Feval, Jeval.

Feval.m - evaluates the force balance equations at a point.

Jeval.m - evaluates the jacobian of the force balance equations at a point.

FBstress.m - This outputs the values of the force balance equations and its jacobian.
Calls:
-Feval.m
-Jeval.m

solve_FBstress.m - solves for the states of the FREE given a specific input pressure by calling lsqnonlin on the function FBstress.

sim_FBstress.m - simulates the steady state deformation of FREE over the interval of input pressures, defined as [params.Pmin, params,Pmax].
Calls:
-solve_FBstress.m

plots_FBstress.m - simulates the steady state deformation of FREE over the interval of input pressures, defined as [params.Pmin, params,Pmax], and at several torsional loads. Then, plots the simulation results over measured data point.
Calls:
-solve_FBstress.m

elastomerModulus.m - symbolically/numerically solves for the elastomer moduli E,G in terms of states of FREE.

elastomerSysID.m - fits lines to elastomer moduli E,G based on measured data points.
Calls:
-elastomerModulus.m




