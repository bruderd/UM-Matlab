function x = setfinalvalues(x, t, y)
% x = setfinalvalues(x, t, y)
%
% Performs the assignment x(...,(end - t + 1):end) = y with the ellipsis being
% the appropriate number of colons (as would be very easy in numpy but not in 
% Matlab because in Matlab, N-dimensional arrays are second-class citizens).
xsize = size(x);
tkeep = xsize(end) - t;
istart = prod(xsize(1:end-1))*tkeep + 1;
xsize(end) = t;
if ~isscalar(y) && ~isequal(size(y), xsize)
    warning('y is the incorrect size. Assignment will fail!');
end
x(istart:end) = y(:);
end%function

