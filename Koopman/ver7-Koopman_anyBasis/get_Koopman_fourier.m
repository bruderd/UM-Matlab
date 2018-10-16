function U = get_Koopman_fourier( x,y, params )
%get_Koopman: Generate finite dimensional approx. of Koopman operator from
%snapshot pairs
%   note: polyLift must be defined for the appropriate system before
%   calling this function.

[n, p] = deal(params.n, params.p);

Px = zeros(length(x), params.N);
Py = zeros(length(x), params.N);
for i = 1:length(x)
    Px(i,:) = fourierLift( x(i, 1:n)', x(i, n+1:n+p)' )';
    Py(i,:) = fourierLift( y(i, 1:n)', y(i, n+1:n+p)' )';
%     Px(i,:) = sinLift( x(i, 1:n)', x(i, n+1:n+p)' )';
%     Py(i,:) = sinLift( y(i, 1:n)', y(i, n+1:n+p)' )';
end

% Better condition Px
c = 1e-3;
% Px = Px + c*eye(size(Px));
Pxinv = dinv(Px, c, params);
% Pxinv = pinv(Px);

disp(cond(Px)); % for RAM debug
U = Pxinv * Py;    % switched to \ to hopefully make computation faster

end