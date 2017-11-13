function q = x2q( x, params )
%x2q: Converts module state to actuator states
%   Detailed explanation goes here

p = params.p;   % number of modules
n = params.n;   % number of actuators in each module (a vector)

q = zeros(2*sum(n),1);
for i = 1:p
    xi = x(1+6*(i-1) : 6*i, 1);     % state of ith module
    [psi, theta, phi] = deal(xi(4), xi(5), xi(6));      % orientation of ith module
    L = params.L(i);        % length of ith module
    
    for j = 1:n(i)
        free_index = sum(n(1:i-1)) + j;     % index of the jth actuator of the ith module
        [a, b] = deal(params.attach(free_index,1), params.attach(free_index,2));   % coordinates of the attachment point
        
        s = -L + sqrt( (L - sin(theta)*a + cos(theta)*sin(psi)*b)^2 + (a^2 + b^2)*phi^2 );
        w = phi;
        
        q(2*free_index - 1 : 2*free_index) = [s, w]';
%         q(2*sum(n(1:i-1)) + 2*j-1:2*j) = [s, w]';
    end
end


end

