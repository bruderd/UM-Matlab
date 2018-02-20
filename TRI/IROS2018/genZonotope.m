function [ zntp, vx, vy ] = genZonotope( X )
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
vy = V(:,6);

zntp = convhull(vx,vy);

% % Plot the zonotope
% figure
% hold on
% plot(X(1,:),X(2,:),'r*');
% plot(vx(zntp), vy(zntp), 'b-');
% hold off

end

