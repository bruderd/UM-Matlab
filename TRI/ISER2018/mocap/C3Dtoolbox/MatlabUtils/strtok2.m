function words = strtok2(str, delim)

words = {};
posDelim = findstr(str, delim);
posDelim = [posDelim length(str)+1];
l = length(posDelim);
currPos = 1;

for i = 1:l
    words{i} = str(currPos:posDelim(i)-1);
    currPos = posDelim(i)+1;
end
