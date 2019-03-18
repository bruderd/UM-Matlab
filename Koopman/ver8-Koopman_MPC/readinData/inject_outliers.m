function [ fucked_data ] = inject_outliers( savepath )
%inject_outliers: Introduce outliers to a bunch of sets of data
%   Detailed explanation goes here


[file,path] = uigetfile('MultiSelect','on');

numFiles = size(file,2);

for i = 1:numFiles
    data = load([ path , file{i}] ); % load the ith data file
    tsteps = size(data.t,1);   % number of timesteps in current data
    maxabs = max( abs( data.x ) );  % maximum absolute value of state

    index = randi( [ 1 , tsteps ] , [randi(floor(tsteps/20)) , 1] );   % choose random indices to impose outliers (no more than 1/20 of points should be corrupted)


    for j = 1:size(index,1)
        data.x(index(j),:) = 2.2*(rand - 0.5) .* maxabs;
    end

    fucked_data{i} = data;
end

% save the new data files if location is provided
if exist('savepath' , 'var')
    for k = 1 : numFiles
        datanow = fucked_data{k};
        x = datanow.x;
        u = datanow.u;
        t = datanow.t;
        newname = [ erase( file{k} , '.mat' ) , '_wOL.mat' ];
        save([ savepath , filesep , newname ] , 't' , 'x' , 'u');
    end 
end


end

