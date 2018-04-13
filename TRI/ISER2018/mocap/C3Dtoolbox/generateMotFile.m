% Generate a motion *.mot file readable by OpenSim
% 
% Tim Dorn
% November 2008
% 
% --------------------------------------------------------------------
% Usage: generateMotFile(dataMatrix, colnames, filename)
% --------------------------------------------------------------------
% 
% Inputs:   dataMatrix = data matrix to write to file
%                       (first column should be time)
%           colnames = cell array of column name strings
%           filename = string containing the output filename (must include extension)
% 
% Outputs:  output motion file
% 
% 
% Notes:    Number of data columns must match the number of column names or
%           an exception will be thrown.
% 
% ----------------------------------------------------------------------

function generateMotFile(dataMatrix, colnames, filename)

[datarows, datacols] = size(dataMatrix);
time = dataMatrix(:,1);
range = [time(1), time(end)];

if length(colnames) ~= datacols
    error('Number of column names do not match the number of columns\n');
end


% MOT File Header
% ---------------

fid = fopen(filename, 'w');
if fid < 0
    fprintf('\nERROR: %s could not be opened for writing...\n\n', filename);
    return
end

fprintf(fid, '%s\nnRows=%d\nnColumns=%d\n\n', filename, datarows, datacols);
fprintf(fid, 'name %s\ndatacolumns %d\ndatarows %d\nrange %f %f\nendheader\n', ...
    filename, datacols, datarows, range(1), range(2));


% MOT File Body
% -------------
cols = [];
for i = 1:datacols,
    if i == 1
        cols = [cols, colnames{i}];
    else
        cols = [cols, sprintf('\t%s', colnames{i})]; 
    end
end
cols = [cols, '\n'];
fprintf(fid, cols);

for i = 1:datarows,
    fprintf(fid, '%20.10f\t', dataMatrix(i,:));
    fprintf(fid, '\n');
end

fclose(fid);
fprintf('Saved motion file: %s\n', filename);

