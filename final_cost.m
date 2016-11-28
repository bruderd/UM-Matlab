function [fc, dfcdx, dfcdu] = final_cost(x, u)

    fc = -10*(sum(x)*sum(u));
    dfcdx = -[10*sum(u) 10+sum(u)];
    dfcdu = -[10*sum(x)];

end