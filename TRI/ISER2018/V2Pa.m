function pressure_Pa = V2Pa( pressure_V )
%V2Pa - converts pressure measurements from the monitor signal of the TR
%pressure regulator to Pa
%   Note: [1 10] V = [0 150] psi

pressure_Pa = (pressure_V * 15) * 6894.76;

end