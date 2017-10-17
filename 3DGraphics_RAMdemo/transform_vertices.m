function vTrans = transform_vertices(v,dirCosine,translation)
% function vTrans = transform_vertices(v,dirCosine,translation)
%
% This function transforms the coordinates of the vertices given in 'v'.
% 'dirCosine' is a rotation 3 x 3 matrix, 'translation' is a translational
% 3-vector. Both are applied to every element in 'v'.
% 'v' is a matrix containing vertices, as they are used in patch objects.
% The reutrn value vTrans contains the coordinates of the transformed
% vertices.
%
%   Nov. 12, 2006 (David Remy, remy@wisc.edu)
%
%   MATLAB Version 7.1

if isempty(v)
        vTrans = [];
        return
    end
    % rotation
    vTrans = (dirCosine*v')';
    % translation
    vTrans = vTrans + repmat(translation,size(vTrans,1),1);
end
