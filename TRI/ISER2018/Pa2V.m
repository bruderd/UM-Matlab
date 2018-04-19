function pressure_V = Pa2V( pressure_Pa, params )
%Pa2V - converts pressure in Pascals to pressure in volts based on the
%current TR settings (which are defined in setParams)
%   Detailed explanation goes here

TRpsimax = params.TRpsimax;             % max pressure [psi]
TRpamax = (1/0.145038) * 1e3 * TRpsimax;     % max pressure [Pa]

pressure_V = min( pressure_Pa * (10/TRpamax), 10 ); % makes sure nothing exceeds 10V

end

