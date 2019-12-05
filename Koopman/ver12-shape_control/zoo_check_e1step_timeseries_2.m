% zoo_check_e1step_timeseries_2

function [ e_calc , e_real ] = zoo_check_e1step_timeseries_2( zx_t , ydes_t , e_t , sys )

for i = 1 : size( ydes_t , 1 ) - ( 1 + sys.params.nd )
    
    e_calc(i,:) = zoo_check_e1step_2( zx_t(i,:)' , ydes_t(i+1,:)' , sys );    % requires same input in both arguments
    e_real(i,:) = e_t(i,:);
    
end

end