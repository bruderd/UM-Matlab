% Find positions of NaN in a 2D vector
% Tim Dorn
% 6th Feb, 2007
% 
% --------------------------------------------------------------------
% Usage: [x,y] = findNaN(var)
% --------------------------------------------------------------------
% 
% Input:   var = variable to find the nans in 
% Output:  returns the (x,y) position of the NaN
% 
% --------------------------------------------------------------------

function [x,y] = findNaN(var)

x = []; y = [];

[m,n] = size(var);

for i = 1:m
    for j = 1:n
        if isnan(var(i,j)),
            x = [x ; j];
            y = [y ; i];
%             disp(['NaN found at position (', num2str(j), ',', ...
%                 num2str(i), ')'])
        end
    end
end

