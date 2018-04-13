% function A = getIndices(cellMain, cellSubset)
% Returns indices of cells containing strings

function A = getIndices(cellMain, cellSubset)

l1 = length(cellSubset);
l2 = length(cellMain);
% A = -1 * ones(1,l1);
A = [];

for i = 1:l1
    for j = 1:l2
        if strcmpi(cellSubset{i}, cellMain{j})
            A = [A j];
%             A(i) = j;
        end
   end
end

return

