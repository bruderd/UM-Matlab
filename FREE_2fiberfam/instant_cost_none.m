function [ic, dicdx, dicdu] = instant_cost_none(x, u, params)

    ic = 0;
    dicdx = zeros(1, length(x));   %must make dimensions consistent, not just 1-D. Although for this special case it should still work
    dicdu = zeros(1, length(u));

end