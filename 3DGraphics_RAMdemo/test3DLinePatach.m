close all
figure
hold on
load('GallopBranch.mat');
X = GallopBranch;
X = diag([1/10,1,5])*X;

%% Just use plot
plot3(X(1,:),X(2,:),X(3,:));
grid on
box on
axis equal
view(3)

%% do a 3D representation
[f, v, n] = create3DLinePatch(X,0.01);
linePatch = patch('faces', f, 'vertices' ,v,'VertexNormals',n);

% Ensure equal scaling in all directions:
daspect([1 1 1])

% Set the face color, remove edges:
set(linePatch, 'FaceColor',[1 0 0]); 
set(linePatch, 'EdgeColor','none');

% Add a default light, right above the current camera position
camlight right 
% Pick renderer
set(linePatch, 'FaceLighting','phong');               



