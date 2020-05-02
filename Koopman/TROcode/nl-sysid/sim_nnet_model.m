function y = sim_nnet_model( u , y0 , u0del , y0del )
%sim_nnet_model
%   Simulates the response of the neural network model (nnet_model) 
%   under the sequence of inputs u.
%   
%   Inputs:
%       u = Tx3 matrix, where T is number of timesteps in trial
%       y0 = 1x2 initial output
%       y0del = 1x2 value of output 1 timestep before start of trial
%       u0 = 1x3 value of input 1 timestep before start of trial

y = zeros( size(u,1) , size(y0,2) );    % preallocate
y(1,:) = y0;
[ y(2,:) , ~ , ~] = nnet_model( u(1,:)' , y(1,:)' , u0del' , y0del' );
for i = 2 : size(u,1)
    [ y(i+1,:) , ~ , ~ ] = nnet_model( u(i,:)' , y(i,:)' , u(i-1,:)' , y(i-1,:)' );
end

end