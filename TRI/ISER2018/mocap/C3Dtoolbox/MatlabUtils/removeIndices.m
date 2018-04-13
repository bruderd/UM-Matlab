function new = removeIndices(cellArray, ind)

n = 1;
new = {};
l = length(cellArray);
for i = 1:l
    if ~any(ind==i)
        new{n} = cellArray{i};
        n = n+1;
    end
end

