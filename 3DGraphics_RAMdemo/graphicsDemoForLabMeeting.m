%
%   Example script to illustrate Matlab 3D graphics
%
%
%   Nov. 12, 2006 (David Remy, remy@wisc.edu)
%
%   MATLAB Version 7.1


fig = clf;
axis off;
set(fig,'Color',[1,1,1]);

% load...
% - faces (triangles), 
% - vertices (3D points), and
% - normals (at these points) of a SIMM-bone file:
[f,v,n] = load_boneFile('foot.asc');
% create a patch object with these faces, vertices and normals:
bonePatch = patch('faces', f, 'vertices' ,v,'VertexNormals',n);

% set up a camera
campos([1,1,1]);        % Define where the camera is   
camtarget([0,0,0]);     % Define where to look at
camup([0,1,0]);         % Define the upwards direction
camproj('perspective'); % Define the projection method
camva(10);              % Define the camera opening angle

% Ensure equal scaling in all directions:
daspect([1 1 1])


% Set the face color:
set(bonePatch, 'FaceColor',[0.90 0.75 0.50]); 

% Add a default light, right above the current camera position
camlight right 

% Choose a good (and slow) renderer
set(bonePatch, 'FaceLighting','gouraud');               
% Choose an even better (and slower) renderer
set(bonePatch, 'FaceLighting','phong');               

% Set the edge color to none, so we don't see them
set(bonePatch, 'EdgeColor','none');

% create the faces and vertices for a cosys representation:
[f,v] = create_cosys_patch(0.1);
% create a patch object with these vertices...
cosysPatch = patch('faces', f, 'vertices' ,v);
set(cosysPatch, 'FaceColor',[1,0,0]); 
set(cosysPatch, 'EdgeColor','none');
set(cosysPatch, 'FaceLighting','phong');               
set(cosysPatch, 'SpecularColorReflectance',[0]);

% let's play a bit with the different material properties:

% set all lights to 0, zo we can see how they work:
set(bonePatch,'SpecularStrength',0,'AmbientStrength',0,'DiffuseStrength',0)

% Ambient light is equal from all directions. Increasing/Decreasing doesn't
% affect the depth impression, the image just get's brighter...
set(bonePatch,'AmbientStrength',1)
set(bonePatch,'AmbientStrength',0)
% Diffuse light depends only on the relative position of the light source
% and the polygon. It is lit most, when the light is coming directly from
% above, less if the light hits at an angle, and it's dark if the light is
% coming from the back. This adds a lot of depth information to the image.
set(bonePatch,'DiffuseStrength',1)
set(bonePatch,'DiffuseStrength',0)
% Specular lights depends on both, the viewer and the light source. It is
% basically the reflection of the light source on the object.
% specular light gives the impression of a very shiny object
set(bonePatch,'SpecularStrength',1)
% The color of the glossy spot depends on the material. 
% It's normally the color of the light source... 
set(bonePatch,'SpecularColorReflectance',[1]);
% .. only metals add their own color
set(bonePatch,'SpecularColorReflectance',[0]);

% The size of the gloosy spot is determined by the Specular exponant:
% The smaller the exponent... the bigger the spots
set(bonePatch,'SpecularExponent',1)
% and vice versa:
set(bonePatch,'SpecularExponent',20)

% Let's set all these values back to a set that looks nice
set(bonePatch,'AmbientStrength',0.4)
set(bonePatch,'DiffuseStrength',0.6)
set(bonePatch,'SpecularStrength',0.2)
set(bonePatch,'SpecularColorReflectance',[0.8]);
set(bonePatch,'SpecularExponent',10)

% another nice trick is to make stuff transparent
set(bonePatch,'FaceAlpha',0.7)
	
% the current properties are listed with the command
get(bonePatch)

% calling set with no property, shows a list of all properties and the
% possible values for them...
set(bonePatch)

% This example is the visualization of a
% Matlab function, with the methods for patch objects:

% create a simple 3D function, and save the X,Y and Z values:
[x,y,z] = peaks();
% Transform this into a patch-object. The last argument is color, for which
% we use the z-value (that gets transformed into a color using a colormap)
[f,v,c] = surf2patch(x,y,z,z);
% To scale this object down to the size of the foot, and to switch the Z
% and the Y axis, we transform every single vertex:
v = transform_vertices(v,[0.03,0,0;0,0,0.01;0,0.03,0],[0,0,0]);
% display as before, the only difference is now, that we also have a
% color value for every vertex:
surfPatch = patch('faces', f, 'vertices' ,v,'CData',c);

% set the other properties
set(surfPatch, 'FaceColor','flat');
set(surfPatch, 'FaceColor','interp');
set(surfPatch, 'EdgeColor','interp');
set(surfPatch, 'FaceLighting','phong');     
set(surfPatch,'AmbientStrength',0.3)
set(surfPatch,'DiffuseStrength',0.6)
set(surfPatch,'SpecularStrength',0.3)
set(surfPatch,'SpecularColorReflectance',0.8);
set(surfPatch,'SpecularExponent',10)

% Finally animate the figure, and save the animation into a .avi video

% set up the output movie:
for r = 0:0.1:2*pi % for a couple of frames
    % the y-direction is periodically scaled:
    scaleY = sin(r);
    % transform the vertices with the oscillating scaling factor:
    v1 = transform_vertices(v,[1,0,0;0,scaleY,0;0,0,1],[0,0,0]);
    % set the vertices of the patch object to the newly calculated ones:
    set(surfPatch, 'vertices',v1);
    % enforce a redraw.
    drawnow;
    
    % and save the output on the screen to a avi video:
    set(fig,'Position',[100,100,800,700]); % to make sure all movieframes are of the same size
end
