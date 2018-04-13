% BRIEF STATMENT ABOUT M FILE
% Tim Dorn
% 29th Jan, 2007
% 
% --------------------------------------------------------------------
% Usage: OUT = FUNCTION(IN1, IN2)
% --------------------------------------------------------------------
% 
% Inputs:   IN1 = 
%           IN2 = 
% 
% Outputs:  OUT = 
% 
% Notes
% -----
% 
% 
% 
% --------------------------------------------------------------------

function OUT = FUNCTION(IN1, IN2)

usage = 'Usage: OUT = FUNCTION(IN1, IN2)';

% ------------------
% ADDITIONAL OPTIONS
% ------------------

% close all

lw = 2;                 % Light plot line width
titlsize = 13;          % Title font size
plottype = 0;           % 0: use subplot windows (for computer viewing)
                        % 1: use individual windows (for image exporting)
saveImages = 1;         % 0: plot only to screen
                        % 1: plot to screen and file 

                        
% Set up some initial parameters and do some initial checks
% ---------------------------------------------------------

if nargin ~= 2,
    disp(usage)
    return
end