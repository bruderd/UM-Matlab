function u = manipulator_input( t , params )
%manipulator_input: input torque into the double pendulum
%   Detailed explanation goes here

% ramp between random values every steplen
u = zeros(params.p , 1);
for i = 1:params.p
    tp = t + (i/params.p) * params.steplen;
    ulast = params.steps( floor(tp/params.steplen) + 1, i )';
    unext = params.steps( ceil(tp/params.steplen) + 1, i )';
    u(i) = ulast + (unext - ulast) .* ( mod(tp,params.steplen) / params.steplen );
end

end

