% iddata2data4sysid
% Converts iddata from a sid workspace to a data4sysid object which can be
% used to identify Koopman models with later versions of sysid code I've
% written.
%
% It doesn't save anything to file, just creates a workspace variable,
% saving will be up to the user

%% training data

% USER NEEDS TO IMPORT THE IDDATA OBJECT CONTAINING ALL THE TRAINING DATA
%   INTO THE WORKSPACE AND CALL IT data4sysid_merged

Ts = data4sysid_merged.Ts{1};   % sampling time

train = cell( 1 , length( data4sysid_merged.y ) - 5 );
for i = 1 : length( data4sysid_merged.y ) - 5
    train{i}.y = data4sysid_merged.y{i}(:,1:3); % don't want derivatives
    train{i}.u = data4sysid_merged.u{i};
    train{i}.t = ( 0 : Ts : ( Ts * (size( train{i}.y , 1 ) - 1) ) )';
end

%% validation data

% USER NEEDS TO IMPORT THE IDDATA OBJECT CONTAINING ALL THE VALIDATION DATA
%   INTO THE WORKSPACE AND CALL IT data4val_merged

val = cell( 1 , length( data4val_merged.y ) );
for i = 1 : length( data4val_merged.y )
    val{i}.y = data4val_merged.y{i}(:,1:3); % don't want derivatives
    val{i}.u = data4val_merged.u{i};
    val{i}.t = ( 0 : Ts : ( Ts * (size( val{i}.y , 1 ) - 1) ) )';
end

%% put em' together!

data4sysid.train = train;
data4sysid.val = val;