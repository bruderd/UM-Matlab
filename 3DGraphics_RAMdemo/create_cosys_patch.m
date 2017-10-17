function [f,v] = getCoSysPatch(cosysSize)
% function [f,v] = getCoSysPatch(cosysSize)
%
% Returns the faces and vertices that form a patch visulization of a
% coordinate system.
%
% The cosys will be drawn in the origin, and aligned with the x,y, and z
% axis
%
% Every axis is cosysSize long.
%
%   Nov. 12, 2006 (David Remy, remy@wisc.edu)
%
%   MATLAB Version 7.1



    % Base sphere at origin
    [x y z] = sphere; 
    [f,v,c] = surf2patch(x,y,z,z);
    v = v*0.05;

    % X Axis
    [x,y,z] = cylinder([0.02,0.02]);
    [f1,v1,c] = surf2patch(x,y,z,z);
    [v,f]=add_patches(v,f,v1,f1);

    [x,y,z] = cylinder([0.4,0]);
    [f1,v1,c] = surf2patch(x,y,z,z);
    v1=v1*0.1;
    v1 = v1 + repmat([0,0,1],size(v1,1),1);
    [v,f]=add_patches(v,f,v1,f1);

    % Y Axis
    [x,y,z] = cylinder([0.02,0.02]);
    [f1,v1,c] = surf2patch(x,y,z,z);
    v1=([1,0,0;0,0,1;0,1,0]*v1')';
    [v,f]=add_patches(v,f,v1,f1);

    [x,y,z] = cylinder([0.4,0]);
    [f1,v1,c] = surf2patch(x,y,z,z);
    v1=v1*0.1;
    v1 = v1 + repmat([0,0,1],size(v1,1),1);
    v1=([1,0,0;0,0,1;0,1,0]*v1')';
    [v,f]=add_patches(v,f,v1,f1);

    % Z Axis
    [x,y,z] = cylinder([0.02,0.02]);
    [f1,v1,c] = surf2patch(x,y,z,z);
    v1=([0,0,1;0,1,0;1,0,0]*v1')';
    [v,f]=add_patches(v,f,v1,f1);

    [x,y,z] = cylinder([0.4,0]);
    [f1,v1,c] = surf2patch(x,y,z,z);
    v1=v1*0.1;
    v1 = v1 + repmat([0,0,1],size(v1,1),1);
    v1=([0,0,1;0,1,0;1,0,0]*v1')';
    [v,f]=add_patches(v,f,v1,f1);
    
    v = transform_vertices(v,eye(3)*cosysSize,[0,0,0]);
    
end