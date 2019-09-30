% zoo_check_pseudoinput_nu

% Ne = size( sysid_unl.traindata.nu , 1 );
Ne = size( sysid_unl.valdata{1}.nu , 1 );
n = size( sysid_unl.traindata.nu , 2 );
m = size( sysid_unl.traindata.u , 2 );
nd = mpc.params.nd;

upseudo_mpc = zeros( Ne-1 , n );
upseudo_real = zeros( Ne-1 , n );
upseudo_diff = zeros( Ne-1 , n );
u_mpc = zeros( Ne-1 , m );
u_real = zeros( Ne-1 , m );
u_diff = zeros( Ne-1 , m );
u_fromnnet = zeros( Ne-1 , m );

% Using traindata
for i = 1 : Ne-nd
    now = i+nd;
    current.y = sysid_unl.traindata.y(i:i+nd,:);
    current.u = sysid_unl.traindata.u(i:i+nd,:);
    current.unl = sysid_unl.traindata.nu(i:i+nd,:);
    
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
    upseudo_real(i,:) = sysid_unl.traindata.nu(i+1,:);
    upseudo_diff(i,:) = abs( Unl(2,:) - sysid_unl.traindata.nu(i+1,:) ); 
    
    % compare the real input to the input computed from mpc pseudoinput
    u_mpc(i,:) = U(2,:);
    u_real(i,:) = sysid_unl.traindata.u(i+1,:);
    u_diff(i,:) = abs( u_mpc(i,:) - u_real(i,:) );
    
    % sanity check: make sure neural network is actually working on perfect pseudoinput data
    u_fromnnet(i,:) = mpc.params.NLinput(sysid_unl.traindata.nu(i,:)')';  % convert to regular input
    
end


% % Using valdata
% % upseudo_mpc_0 = sysid_unl.valdata{1}.nu(1:1+nd,:);  % DONT USE PERFECT PSEUDOINPUT
% for i = 1 : Ne-nd
%     now = i+nd;
%     current.y = sysid_unl.valdata{1}.y(i:i+nd,:);
%     current.u = sysid_unl.valdata{1}.u(i:i+nd,:);
%     current.unl = sysid_unl.valdata{1}.nu(i:i+nd,:);
% %     if i == 1   % DONT USE PERFECT PSEUDOINPUT
% %         current.unl = upseudo_mpc_0;
% %     else
% %         current.unl = upseudo_mpc(i-1,:);
% %     end
%     
%     Nh = mpc.horizon;
%     if now + Nh < Ne
%         refhor = sysid_unl.valdata{1}.y( now : now+Nh , : );
%     elseif now < Ne
%         refhor = sysid_unl.valdata{1}.y( now : end , : );
%         refhor = [ refhor ; kron( ones(Ne-now,1) , sysid_unl.valdata{1}.y(end,:) ) ];
%     else
%         refhor = kron( ones(Nh,1) , sysid_unl.valdata{1}.y(end,:) );
%     end
%     
%     % find the pseudoinput that mpc calculates
%     [ Unl , z ] = mpc.get_mpcInput_noIC( current , refhor );
%     U = mpc.params.NLinput(Unl')';  % convert to regular input
%     
%     % compare the pseudoinput to the actual error
%     upseudo_mpc(i,:) = Unl(1,:);
%     upseudo_real(i,:) = sysid_unl.valdata{1}.nu(i+1,:);
%     upseudo_diff(i,:) = abs( Unl(2,:) - sysid_unl.valdata{1}.nu(i+1,:) ); 
%     
%     % compare the real input to the input computed from mpc pseudoinput
%     u_mpc(i,:) = U(2,:);
%     u_real(i,:) = sysid_unl.valdata{1}.u(i+1,:);
%     u_diff(i,:) = abs( u_mpc(i,:) - u_real(i,:) );
%     
%     % sanity check: make sure neural network is actually working on perfect pseudoinput data
%     u_fromnnet(i,:) = mpc.params.NLinput(sysid_unl.valdata{1}.nu(i,:)')';  % convert to regular input
%     
% end

