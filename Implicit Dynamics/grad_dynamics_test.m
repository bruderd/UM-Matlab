function grad_coneq = grad_dynamics_test(s, params)

%     [~, ~, ~, grad_coneq] = dynamics_noinput(s, params);
    
    [~, ~, ~, grad_coneq] = dynamics_winput(s, params);

end