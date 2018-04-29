function pressure_Pa = V2Pa( pressure_V )
%V2Pa - converts pressure monitor signals from TR pressure regulator to Pa
%   Note: [0 10] V = [0 150] psi

pressure_Pa = pressure_V * (150/10) * 6894.76;

end