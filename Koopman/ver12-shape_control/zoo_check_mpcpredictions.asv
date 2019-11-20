function [ Empc , Ereal ] = zoo_check_mpcpredictions( index , sysid_class , mpc_class )
% zoo_check_mpcpredictions
%
% Check the accuracy of the pseudoinput chosen by MPC over a finite time
% horizon by comparing it to training data.
%
% must have 'sysid_unl' and 'mpc' classes defined before calling this.
%

% set initial condition for this step
traj.y = sysid_class.traindata.y( index , : );
traj.u = sysid_class.traindata.u( index , : );
traj.e = sysid_class.traindata.e( index , : );

ref = sysid_class.traindata.y( index : index + mpc_class.horizon , : );
Ereal = sysid_class.traindata.e( index : index + mpc_class.horizon , : );

[ Empc , zx ]= mpc_class.get_mpcInput_unl( traj , ref );

% calculate the best mpc input for single time step
z0 = sysid_class.lift.zx( traj.y' );
zdes = sysid_class.lift.zx( ref(2,:)' );
C = eye( sysid_class.params.nzx );
d = -sysid_class.model.A * z0 + zdes;

enow = lsqlin( C , d , [] , [] );

end