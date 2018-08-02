function u = get_u( t )
%get_u: Defines input into the system as a function of time
%   u is a frequency/amplitude varying sinusoidal input for generating sysid data

% u = sin(2*pi*(2*t-cos(t)));
u = 4*sin( (1/(2*pi)) * t) .* sin( 3*t - 1.5*cos(t) );


end

