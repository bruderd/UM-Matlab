function theta = alpha2theta(alpha)
%alpha2theta: Converts relative joint angles (alpha) to absolute joint
%angles (theta).
%   Detailed explanation goes here

theta = zeros( size(alpha) );

% if input is symbolic, so should output
if isa( alpha , 'sym' )
    theta = sym(theta);
end

for i = 1 : length(alpha)
   theta(i) =  sum(alpha(1:i));
end

end

