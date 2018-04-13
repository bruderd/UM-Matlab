function s = removeSpaces(str)
s = [];
l = length(str);

for i = 1:l
	if ~isspace(str(i))
       s = [s str(i)]; 
    end
end

