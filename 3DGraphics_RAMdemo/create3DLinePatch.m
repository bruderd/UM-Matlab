function [faces, vertices, normals] = create3DLinePatch(X,radius)
% returns the vertices and faces of a three dimensional tubular 'line' with
% the given radius that runs along the line defined by vectors in X, which
% is a 3xm matrix 

    m = size(X,2);
    n_segments = 10; % segments around circumference
    %compute local directions u,v
    u = zeros(3,m);
    v = zeros(3,m);
    % First element:
    d = X(:,2)-X(:,1);  % direction
    n = rand(3,1);      % normal (should remain the same along trace)
    u(:,1) = cross(d,n);
    v(:,1) = cross(d,u(:,1));
    u(:,1) = u(:,1)/norm(u(:,1));
    v(:,1) = v(:,1)/norm(v(:,1));
    % along trace
    for i = 2:m-1
        d = X(:,i+1)-X(:,i-1);  % direction
        n = -v(:,i-1);          % normal (should remain the same along trace)
        u(:,i) = cross(d,n);
        v(:,i) = cross(d,u(:,i));
        u(:,i) = u(:,i)/norm(u(:,i));
        v(:,i) = v(:,i)/norm(v(:,i));
    end 
    % last element:
    d = X(:,m)-X(:,m-1);  % direction
    n = -v(:,m-1);      % normal (should remain the same along trace)
    u(:,m) = cross(d,n);
    v(:,m) = cross(d,u(:,m));
    u(:,m) = u(:,m)/norm(u(:,m));
    v(:,m) = v(:,m)/norm(v(:,m));
    
    % set up vertices and normals:
    vertices = zeros(3,m*n_segments+2);
    normals = zeros(3,m*n_segments+2);
    % Endcaps
    vertices(:,1) = X(:,1);
    normals(:,1) = X(:,1)-X(:,2);
    vertices(:,end) = X(:,end);
    normals(:,end) = X(:,m)-X(:,m-1);
    for i = 1:m
        for j = 1:n_segments
            vertices(:,vertexIndex(i,j)) = X(:,i) + radius*u(:,i)*sin(j/n_segments*2*pi) + radius*v(:,i)*cos(j/n_segments*2*pi);
            normals(:,vertexIndex(i,j))  =                 u(:,i)*sin(j/n_segments*2*pi) +        v(:,i)*cos(j/n_segments*2*pi);
        end
    end 
    
    %set up faces:
    faces = zeros(4,(m+1)*n_segments);
    for j = 1:n_segments
        faces(:,j) = [1;vertexIndex(1,j);vertexIndex(1,j+1);1];
        faces(:,m*n_segments+j) = [m*n_segments;vertexIndex(m,j);vertexIndex(m,j+1);m*n_segments];
    end
    for i = 1:m-1
        for j = 1:n_segments
            faces(:,i*n_segments+j) = [vertexIndex(i,j);vertexIndex(i+1,j);vertexIndex(i+1,j+1);vertexIndex(i,j+1)];
        end
    end 
    
    % transpose everything:
    faces = faces';
    vertices = vertices';
    normals = normals';
    
    % helper function
    function index = vertexIndex(i_,j_)
        if j_ > n_segments
            j_ = j_-n_segments;
        end
        index = (i_-1)*n_segments+j_+1;
    end
end