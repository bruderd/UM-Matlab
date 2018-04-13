function ACHANNEL = getanalogchannel(itf, signalname, index1, index2)
% GETANALOGCHANNEL - returns nx1 array containing analog data.  The
% returned analog data is scaled with offsets removed, but is in the force plate
% coordinate system.
% 
%   USAGE:  ACHANNEL = getanalogchannel(itf, signalname, index1*, index2*)
%           * = not a necessary input
%   INPUTS:
%   itf         = variable name used for COM object
%   signalname  = string name of desired signal
%   index1      = start frame index, all frames if not used as an argument
%   index2      = end frame index, all frames if not used as an argument
%   OUTPUTS:
%   ACHANNEL    = nx1 matrix of analog data

%   C3D directory contains C3DServer activation and wrapper Matlab functions.
%   This function written by:
%   Matthew R. Walker, MSc. <matthewwalker_1@hotmail.com>
%   Michael J. Rainbow, BS. <Michael_Rainbow@brown.edu>
%   Motion Analysis Lab, Shriners Hospitals for Children, Erie, PA, USA
%   Questions and/or comments are most welcome.  
%   Last Updated: April 21, 2006
%   Created in: MATLAB Version 7.0.1.24704 (R14) Service Pack 1
%               O/S: MS Windows XP Version 5.1 (Build 2600: Service Pack 2)
%   
%   Please retain the author names, and give acknowledgement where necessary.  
%   DISCLAIMER: The use of these functions is at your own risk.  
%   The authors do not assume any responsibility related to the use 
%   of this code, and do not guarantee its correctness. 
     

if nargin == 2, 
    index1 = itf.GetVideoFrame(0); % frame start
    index2 = itf.GetVideoFrame(1); % frame end
elseif nargin == 1 || nargin == 3,
    disp('Error: wrong number of inputs.');
    help getanalogchannel;
    return
end

% Not sure what this does but the exception it creates can be discarded
% because the program still works... TD: Nov 2008
% ------------------------------------------------
%  fIndex = itf.GetParameterIndex('ANALOG', 'FORMAT');
%  if fIndex == -1,
%      disp('Sorry, this function does not support your analog data format at this time.');
%      return
%  else
%     format = itf.GetParameterValue(fIndex, 0);
%     if upper(format(1)) == 'S',
%         disp('Error: only supports VICON unsigned analog data at this time');
%     return
%     else
% %         disp('C3D file contains unsigned analog data, retrieving data...');
%     end
%  end

 
 nIndex = itf.GetParameterIndex('ANALOG', 'LABELS');
 nItems = itf.GetParameterLength(nIndex);
 signalname = upper(signalname);
 signal_index = -1;
    
for i = 1 : nItems,
    channel_name = removeSpaces(itf.GetParameterValue(nIndex, i - 1));
    newstring = channel_name(1:min(findstr(channel_name, ' '))-1);
    if strmatch(newstring, [], 'exact'),
        newstring = channel_name;
    end
    
%     if strmatch(upper(newstring), signalname, 'exact') == 1,
    if findstr(upper(newstring), upper(removeSpaces(signalname)))
        signal_index = i - 1;
    end  
end

if signal_index ~= -1,
    ACHANNEL = itf.GetAnalogDataEx(signal_index,index1,index2,'1',0,0,'0');
    ACHANNEL = cell2mat(ACHANNEL); 
else
    fprintf('Warning: Analog Channel [%s] not found.\n', upper(signalname));
    ACHANNEL = [];
    return
end




% ADDED BY TIM 12/6/08
% Since this ACHANNEL value represents the analog values at the index1,
% index2 VIDEO frames and they are inclusive, we must subtract (r-1) frames
% from the end to represent the analog data frrom and including index1, to
% the start of index2 (exclusive of it)

r = double(itf.GetAnalogVideoRatio);
l = length(ACHANNEL);
ACHANNEL = ACHANNEL(1:l-(r-1));

%--------------------------------------------------------------------------