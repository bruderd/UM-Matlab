% Change coordinate system between Vicon / Force Plate / 
% OR any other arbitrary coordinate system
% Tim Dorn
% 16th Jun 2008
% 
% --------------------------------------------------------------------
% Usage: newSysData = coordChange(oldSysData, dirVec)
% --------------------------------------------------------------------
% 
% Inputs:   oldSysData = matrix of a multiple of 3 rows in old coord sys
%                       e.g. [x, y, z] OR [x1, y1, z1, x2, y2, z2]
% 
%           dirVec = transformation vector [x, y, z] form
% 
%
% Outputs:  newSysData = matrix of a multiple of 3 rows in new coord sys
%                       e.g. [x, y, z] OR [x1, y1, z1, x2, y2, z2] 
% 
% 
% Notes
% -----
% 
% newSystem = sign(dirVec) * oldSysData(|dirVec|)
% 
% No time vectors should be in the data set as inputs. 
% No time vectors are returned as output.
% 
% --------------------------------------------------------------------

function newSysData = coordChange(oldSysData, dirVec)

usage = 'Usage: newSysData = coordChange(oldSysData, dirVec)';

if nargin ~= 2,
    disp(usage)
    return
end



% Check input data for correctness
% --------------------------------

if isvector(dirVec) == 0 || length(dirVec) ~= 3,
    error('Incorrect input parameters...\n%s', usage)
end

[m,n] = size(oldSysData);
dirTmp = dirVec;
dirVec = [];

if mod(m,3) == 0,
    r = floor(m/3);
    for i = 1:r
        dirVec = [dirVec, sign(dirTmp).*(abs(dirTmp)+(3*(i-1)))];
    end
else
    error('Input data must have a multiple of 3 rows i.e. (x,y,z) / (x,y,z,x,y,z)...\n%s', usage)
end
    


% Perform the transformation
% --------------------------

for i = 1:m
    newSysData(i,:) = sign(dirVec(i)) * oldSysData(abs(dirVec(i)), :);
end


