% append constant text to end of each cell member

function newC = appendEachCell(c, txt)

l = length(c);
newC = {};

for i = 1:l
    newC{i} = [c{i} txt];
end
