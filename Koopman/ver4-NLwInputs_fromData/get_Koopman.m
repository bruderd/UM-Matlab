function U = get_Koopman( x,y, params )
%get_Koopman: Generate finite dimensional approx. of Koopman operator from
%snapshot pairs
%   note: polyLift must be defined for the appropriate system before
%   calling this function.

[n, p] = deal(params.n, params.p);

%Px = zeros(length(x),)
for i = 1:length(x)
    Px(i,:) = polyLift( x(i, 1:n)', x(i, n+1:n+p)' )';
    Py(i,:) = polyLift( y(i, 1:n)', y(i, n+1:n+p)' )';
end

Pxinv = pinv(Px);
disp(cond(Px)); % for RAM debug
U = Pxinv * Py;

end