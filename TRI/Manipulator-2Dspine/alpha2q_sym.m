function q = alpha2q_sym( alpha, params )
%Converts alpha to q
%   Detailed explanation goes here


s = zeros(params.n,1);
s = sym(s);
s(1,1) = params.dl - params.attach(1) * ( (alpha(1) / 2 + alpha(1)^3 / 24 ) );
s(2,1) = params.dl - params.attach(2) * ( (alpha(1) / 2 + alpha(1)^3 / 24 ) );
for i = 1:params.n
    for j = 2:params.p
        partSum = params.dl - params.attach(i) * ( (alpha(j-1) / 2 + alpha(j-1)^3 / 24) + (alpha(j) / 2 + alpha(j)^3 / 24 ) );
        s(i,1) = s(i,1) + partSum; 
    end
end
s = s - params.L * ones(size(s));

q = s;

end

