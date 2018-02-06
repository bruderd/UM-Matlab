function xc = infercolloc(r, x, method)
% xc = infercolloc(r, x, [method])
%
% Interpolates the `[Nx, Nt + 1]` array `x` at the collocation points defined in
% `r` to yield an `[Nx, Nc, Nt]` array of values at the collocation points.
%
% Default for `method` is 'linear'. Other options are 'pchip' and 'spline'. See
% `help interp1` for details.
narginchk(2, 3);
if nargin() < 3
    method = 'linear';
end

Nx = size(x, 1);
Nt = size(x, 2) - 1;

r = r(2:end-1); % Get rid of endpoints.
r = r(:)'; % Make sure it's a row vector.
Nc = length(r);

if Nt == 0
    xc = repmat(x, [1, Nc, 1]);
else
    t = 0:Nt;
    tc = kron(1:Nt, r);
    xc = interp1(t, x', tc)';

    xc = reshape(xc, [Nx, Nc, Nt]);
end

end%function

