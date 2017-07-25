function f2tau = calc_f2tau( x, params )
%calc_Jv
%   Defines the Jacobian, Jv, such that fq = Jv * u

n = params.numFREEs;   % the number of FREE actuators in system
a = params.xattach;   % x-cordinate of FREE attachment points (array)
b = params.yattach;   % y-coordinate of FREE attachment points (array)

f2tau = zeros(3, 2*n);      % initialize Jv
for i = 1:n
    f2tau(1,i) = -b(i);   
    f2tau(2,i) = -a(i);
    f2tau(3,n+i) = 1;
end

end

