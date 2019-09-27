% zoo_check_pseudoinput

Ne = size( sysid_unl.traindata.e , 1 );
n = size( sysid_unl.traindata.e , 2 );
m = size( sysid_unl.traindata.u , 2 );
nd = mpc.params.nd;

upseudo_mpc = zeros( Ne-1 , n );
upseudo_real = zeros( Ne-1 , n );
upseudo_diff = zeros( Ne-1 , n );
u_mpc = zeros( Ne-1 , m );
u_real = zeros( Ne-1 , m );
u_diff = zeros( Ne-1 , m );
u_fromnnet = zeros( Ne-1 , m );
for i = 1 : Ne-nd
    now = i+nd;
    current.y = sysid_unl.traindata.y(i:i+nd,:);
    current.u = sysid_unl.traindata.u(i:i+nd,:);
    current.unl = sysid_unl.traindata.e(i+nd,:);
    
    Nh = mpc.horizon;
    if now + Nh < Ne
        refhor = sysid_unl.traindata.y( now : now+Nh , : );
    elseif now < Ne
        refhor = sysid_unl.traindata.y( now : end , : );
        refhor = [ refhor ; kron( ones(Ne-now,1) , sysid_unl.traindata.y(end,:) ) ];
    else
        refhor = kron( ones(Nh,1) , sysid_unl.traindata.y(end,:) );
    end
    
    % find the pseudoinput that mpc calculates
    [ Unl , z ] = mpc.get_mpcInput( current , refhor );
    U = mpc.params.NLinput(Unl')';  % convert to regular input
    
    % compare the pseudoinput to the actual error
    upseudo_mpc(i,:) = Unl(2,:);
    upseudo_real(i,:) = sysid_unl.traindata.e(i+1,:);
    upseudo_diff(i,:) = abs( Unl(2,:) - sysid_unl.traindata.e(i+1,:) ); 
    
    % compare the real input to the input computed from mpc pseudoinput
    u_mpc(i,:) = U(2,:);
    u_real(i,:) = sysid_unl.traindata.u(i+1,:);
    u_diff(i,:) = abs( u_mpc(i,:) - u_real(i,:) );
    
    % sanity check: make sure neural network is actually working on perfect pseudoinput data
    u_fromnnet(i,:) = mpc.params.NLinput(sysid_unl.traindata.e(i,:)')';  % convert to regular input
    
end