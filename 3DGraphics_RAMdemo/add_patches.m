function [v, f] = addPatches(v1,f1,v2,f2)
% function [v, f] = addPatches(v1,f1,v2,f2)
%
% Assuming that v1, f1, v2, and f2 define vertices and faces of two patch
% objects, this function returns the vertices and faces of a new combined
% patch object.
%
%   Nov. 12, 2006 (David Remy, remy@wisc.edu)
%
%   MATLAB Version 7.1


    f2 = f2 + repmat(size(v1,1),size(f2));
    v = [v1; v2];
    f = [f1; f2];
end