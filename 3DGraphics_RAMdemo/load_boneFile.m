function [fout, vout, nout] = load_boneFile(filename)
% function [fout, vout] = load_boneFile(filename)
%
% Reads a SIMM bonefile (*.asc) and returns two matrices containing the
% faces and vertices as they are used by the patch-function.
% 
% input   filename - the full path and filename of the bone file.
%
% output  vout     - A vector of vertices, i.e 3D points, that define the
%                    positions of the edges of the faces.
%         fout     - A vector of faces (i.e. triangles), defined my the
%                    vertices these triangles span.
%
% Example of how to include this into a Matlab program:
% 
% figure
% [f,v] = load_boneFile('foot.asc');
% p = patch('faces', f, 'vertices' ,v);
% set(p, 'FaceColor',[0.90 0.75 0.50]); % Set the face color
% set(p, 'FaceLighting','phong');       % Choose a good (and slow) renderer
% set(p, 'EdgeColor','none');           % Don't show the edges
% light                                 % Add a default light
% daspect([1 1 1])                      % Ensure equal scaling in alldirections
% view(3)                               % Isometric view
%
%   Nov. 12, 2006 (David Remy, remy@wisc.edu)
%
%   MATLAB Version 7.1

fid=fopen(filename, 'r');               % Open the file, assumes a SIMM bone file in asc ASCII format.
if fid == -1 
    error('File could not be opened, check name or path.')
end

nElements = [0;0];  % Number of Vertices and Faces in the input file
vnum=0;             % Vertex number counter.
fnum=0;             % Face number counter

v=[];
f=[];

while feof(fid) == 0                    % test for end of file, if not then do stuff
    tline = fgetl(fid);                 % reads a line of data from file.
    fword = sscanf(tline, '%s ');       % make the line a character string
    if isempty(fword )                  % ignor empty lines
        continue
    end
    if strncmpi(fword, 'NORM_ASCII',10) == 1;     % Checking if it's the first line "NORM_ASCII"
       continue                                   % ...and ignore it  
    end        
    if (nElements(1) == 0)                        % Read the header
        nElements = sscanf(tline,'%d');          
        tline = fgetl(fid);                       % read the following line (the bounding box) 
        continue                                  % ... and ignore it
    end
    if vnum<nElements(1)                          % read the vertices and store them into v
        vnum = vnum + 1; 
        v(:,vnum) = sscanf(tline, '%f');
        continue
    else                                          % done with vertices 
        if fnum<nElements(2)                      % read the faces.
            fnum = fnum + 1; 
            poly = sscanf(tline, '%d');
            for i=1:poly(1)-2                     % break down faces with more then 3 edges into...
                f=[f,[1+poly(2);1+poly(i+2);1+poly(i+3)]]; %... basic polygones and store them into f
            end
            continue
        end
    end
end
fout(:,1) = f(3,:)';          % Transpose the Matrix for the use in patch().
fout(:,2) = f(2,:)';          % And switch the first and the third element, so that matlab
fout(:,3) = f(1,:)';          % will create normal vectors that point outwards

vout = v(1:3,:)';   % Transpose the matrix, and only use the first three elements. 
nout = v(4:6,:)';   % The next six elements are the normal vector.
fclose(fid);
