function u = get_simInput( t, data, params )
%get_simInput: Estimates the input at time t from the data
%   Used to compare simulation of sysid'd system to performance of real
%   system.

for j = 1 : length(data.u)-1
    if t >= data.t(j) && t <= data.t(j+1)
        u = data.u(j,:)';
        break
    end
end

end

