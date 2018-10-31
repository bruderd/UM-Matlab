function util_plotpos(frame, xcoord, ycoord)
%UTIL_PLOTPOS Plot laser coordinates in the image frame provided.
%
%    UTIL_PLOTPOS(FRAME, X, Y) plots the image FRAME provided and 
%    a set of cross hairs at the X and Y coordinates.
%

%    DH 2-16-03
%    Copyright 2001-2005 The MathWorks, Inc.

% Display the image frame
imshow(frame);

% Superimpose cross hairs.
hold on
plot(xcoord, ycoord, 'yo')
plot([1 size(frame, 2)], [ycoord ycoord], 'y-')
plot([xcoord xcoord], [1 size(frame, 1)], 'y-')
hold off;
