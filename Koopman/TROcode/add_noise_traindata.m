% add noise to TRO robot data: softrobot_train-13_val-4.mat data

% load in the system identifcation data struct
load('C:\Users\danie\OneDrive\Documents\MATLAB\UM-Matlab\Koopman\TROcode\datafiles\softrobot_train-13_val-4.mat')

for i = 1 : length( train )
    % add uniformly distributed noise
%     offset = rand;
%     train{i}.y = train{i}.y + ( 2*randn( size(train{i}.y ) ) - ( 0.5 + offset ) );

    % add normally distributed noise
%     offset = randn( size( train{i}.y ) ); % different offset for each dimension
%     train{i}.y = train{i}.y + ( randn( size(train{i}.y ) ) + offset );
    
    % add another signal instead of random noise
    train{i}.y = train{i}.y + 4*sin( ( 1 : length( train{i}.t ) ) * 1e0 )';
end

% save the noisier data
save( [ 'datafiles' , filesep , 'sofrobot_noisy.mat' ] , 'train' , 'val' );