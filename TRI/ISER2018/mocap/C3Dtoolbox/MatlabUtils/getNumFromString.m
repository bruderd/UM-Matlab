% Extract a number from a string
% Tim Dorn
% 7th Apr 2008
% 
% example: getNumFromString('hello123') = 123
%          getNumFromString('a1b2c3')   = 123
% 
% ---------------------------------------------------

function num = getNumFromString(string)

strNum = [];
l = length(string);

for i = 1:l
    tmp = str2double(string(i));
    if isnan(tmp) || strcmp('i', string(i)) || strcmp('j', string(i))
        % not a valid number
    else
        % is a valid number
        strNum = [strNum, num2str(tmp)];
    end
end

num = str2double(strNum);
