% zoo_check_nu1step_timeseries_2

function [ nu_calc , nu_real ] = zoo_check_nu1step_timeseries_2( zx_t , ydes_t , nu_t , sys )

for i = 1 : size( ydes_t , 1 ) - ( 1 + sys.params.nd )
    
    nu_calc(i,:) = zoo_check_nu1step_2( zx_t(i,:)' , ydes_t(i+1,:)' , sys );    % requires same input in both arguments
    nu_real(i,:) = nu_t(i,:);
    
end

end