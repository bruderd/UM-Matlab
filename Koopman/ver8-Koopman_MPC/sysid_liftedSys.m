function out = sysid_liftedSys( U, params )
%sysid_liftedSys: Defines the A, B, and C matrices that describe the linear
%lifted system z+ = Az + Bu, x = Cz.
%   Detailed explanation goes here

UT = U';    % transpose of koopman operator

A = UT( 1 : params.N , 1 : params.N );
B = UT( 1 : params.N , params.N+1 : end );

if strcmp( params.basisID, 'poly' )
    C = [zeros(params.n , 1) , eye(params.n) , zeros( params.n , params.N - params.n - 1 )];   % if poly we want to skip over first element of lifted state which is "1"
else
    C = [eye(params.ny), zeros(params.ny , params.N - params.ny)];   % C selects the first ny entries of the lifted state (so output can be different than state)
end

%% define outputs
out.A = A;
out.B = B;
out.C = C;
out.sys = ss(A,B,C,0, params.Ts);  % discrete state space system object
out.params = params;    % save system parameters as part of system struct

end

