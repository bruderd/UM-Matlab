% zoo_check_e1step_timeseries

function [ e_calc , e_real ] = zoo_check_e1step_timeseries( zeta_t , ydes_t , e_t , sys )

for i = 1 : size( ydes_t , 1 ) - ( 1 + sys.params.nd )
    e_calc(i,:) = zoo_check_e1step( zeta_t(i,:)' , ydes_t(i+1,:)' , sys );
    e_real(i,:) = e_t(i,:);
end

end