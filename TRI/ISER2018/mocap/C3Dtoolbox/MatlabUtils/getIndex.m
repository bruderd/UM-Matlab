% function i = getIndex(cell, str)
% Returns index of first exact matching string

function i = getIndex(cell, str)

l = length(cell);
for i = 1:l
   if strcmp(str, cell{i})
       return;
   end
end
i = 0;
% fprintf('Could not find %s in cell array\n', str);