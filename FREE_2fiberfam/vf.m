function [f, dfdx, dfdu, df_ddxdt] = vf(x, u, dxdt, params)
% This function selects the proper version of the dynamics to use
%

%     [f, dfdx, dfdu, df_ddxdt] = vf2(x, u, dxdt, params);  
    
    [f, dfdx, dfdu, df_ddxdt] = vf_CASES_v6(x, u, dxdt, params);
       
end
    