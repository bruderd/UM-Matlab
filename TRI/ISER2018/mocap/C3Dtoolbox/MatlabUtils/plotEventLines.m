% Plot vertical event lines
% -------------------------
% all = 1: plot all events (including general events)
% all = 0: plot only foot events (all except general events)
% 
% -------------------------

function plotEventLines(C3Dkey, timeVectorToUse, all)

if nargin == 2
    all = 1;
elseif nargin ~= 3
    disp('Error plotting event lines: USAGE: plotEventLines(C3Dkey, timeVectorToUse, all*\n')
    return
end

col = 'm:';
ratioTxtAbove = 1.05;
timeVec = C3Dkey.event.(timeVectorToUse);

% Only plot events within the time frame of the x axes plot window
xl = xlim;
events2use = [];
for i = 1:length(timeVec)
    if timeVec(i) >= xl(1) && timeVec(i) <= xl(end)
       events2use = [events2use i]; 
    end
end
% disp(events2use)


% Plot event lines
try
    for i = events2use
        % all = 0
        if ~all
            if ~strcmp(C3Dkey.event.txt{i}, 'rGEN') && ...
               ~strcmp(C3Dkey.event.txt{i}, 'lGEN') && ...
               ~strcmp(C3Dkey.event.txt{i}, 'GEN')
                vline(timeVec(i), col, C3Dkey.event.txt{i}, ratioTxtAbove);
            end
            
        else
            vline(timeVec(i), col, C3Dkey.event.txt{i}, ratioTxtAbove);
        end
        
    end
catch
    return
end
