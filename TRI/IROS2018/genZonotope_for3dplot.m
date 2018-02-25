function [ zntp, vx, vy ] = genZonotope_for3dplot( X )
%genZonotope: Constructs zonotope from input vectors
%   Computes convex hull of all binary combinations of vectors in X.
%   X should be a horizontal concatenation of column vectors.

nx = size(X, 1);
N = size(X, 2);

for i = 1:N^2
    m = dec2bin(i-1,N);
    for j = 1:N
        M(i,j) = str2double(m(j));
    end
end

V = M*X';
vx = V(:,3);
vy = V(:,2);
vz = V(:,1);

zntp = convhull(vx,vy,vz);

% Plot the zonotope
figure
hold on
plot3(X(1,:),X(2,:),X(3,:),'r*');
plot3(vx(zntp), vy(zntp), vx(zntp), 'b-');
hold off

% SOMETHING IS WRONG DONT TRUST THIS

end