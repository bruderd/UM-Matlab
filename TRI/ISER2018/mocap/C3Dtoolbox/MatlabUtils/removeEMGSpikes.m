% Remove spikes in EMG and replace with normative data
% Tim Dorn
% 2nd Oct 2008
% 
% --------------------------------------------------------------------
% Usage: out = removeEMGSpikes(inputData)
% --------------------------------------------------------------------
% 
% Inputs:   inputData = input EMG data
% 
% Outputs:  out = output EMG data
% 
% Notes:    Manual intervention is used to remove spikes
% 
% --------------------------------------------------------------------

function out = removeEMGSpikes(inputData)
figure(999)
set(gcf, 'Position', [12 704, 1896, 420]);
plot(inputData)
axis tight


while(1)
    title('Select 2 Source Points for Replacement (ENTER to end)')
    [xsrc, ysrc, but] = ginput(2);

    xsrc = sort(xsrc);

    if isempty(xsrc) || length(xsrc) == 1,
        out = inputData;
        close(999)
        return
    end

    sourceFrames = floor(xsrc(1)):ceil(xsrc(2));
    ShadePlotForEmpahsisVert(xsrc, 'y', 0.3);

    title('Select 2 Destination Points to replace shaded yellow area (ENTER to end)')
    [xdst, ydst, but] = ginput(2);
    xdst = sort(xdst);

    if isempty(xdst) || length(xdst) == 1,
        out = inputData;
        close(999)
        return
    end
    
    destFrames = floor(xdst(1)):ceil(xdst(2));
    tmp = inputData;

    inputData(sourceFrames) = resamp2(inputData(destFrames)', length(sourceFrames));
    hold off
    plot(inputData)
    axis tight
    hold on
    plot(sourceFrames, tmp(sourceFrames), 'r--')
end

