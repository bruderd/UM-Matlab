function pressure_Pa = Vin2Pa( pressure_V, params )
%V2Pa - converts pressure input signals to TR pressure regulator to Pa

pressure_Pa = pressure_V * (params.TRpsimax/10) * 6894.76;

end