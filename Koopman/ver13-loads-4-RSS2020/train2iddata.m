% train2iddata
%   Converts the training/validation data into the iddata format so that it
%   can be used with the system identification toolbox.


% load in the training data
data4sysid = load([ 'datafiles' , filesep , 'softarm_0-300g_3marks_ramp-0p5-3s_trim_train-49_val-49_2020-01-16_11-23.mat' ]);

% Assume timesetp is 12 hz
Ts = 0.083;

% % get average timestep from data
% Ts_i = zeros( size( data4sysd.train ) ); 
% for i = 1 : length( data4sysid.train )
%     Ts_i(i) = mean( data4sysid.train{i}.t(2:end) - data4sysid.train{i}.t(1:end-1) )
% end
% Ts = mean( Ts_i );

% convert all the trainin data to iddata
for i = 1 : length( data4sysid.train )
    train_iddata{i} = iddata( data4sysid.train{i}.y , data4sysid.train{i}.u , Ts);
    val_iddata{i} = iddata( data4sysid.val{i}.y , data4sysid.val{i}.u , Ts);
end

% merge the data trials
train_merged = train_iddata{1};
val_merged = val_iddata{1};
for i = 2 : length( data4sysid.train ) 
    train_merged = merge( train_merged , train_iddata{i});
    val_merged = merge( val_merged , val_iddata{i});
end