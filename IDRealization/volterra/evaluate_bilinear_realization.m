function [ yout, ybl ] = evaluate_bilinear_realization(A,N,B,C,kfinal)


%% simulate real system

% kfinal = 10;
numtrials = 1;

step_magnitude = 5e-1;
utrain = zeros( numtrials , kfinal );
for j = 1 : numtrials
    % generate input as a random walk in [-1,1]
    steps = randi([-1 1],1,kfinal);
    scales = rand(1,kfinal);
    u = zeros( 1 , kfinal);
    u(1) = 1e-2 * rand;
    for i = 2 : length(steps)
        if u(i-1) >= 1
            u(i) = u(i-1) - step_magnitude;
        elseif u(i-1) <= -1
            u(i) = u(i-1) + step_magnitude;
        else
            u(i) = u(i-1) + 2*step_magnitude * scales(i) * steps(i);
        end
    end
    utrain(j,:) = u;
end

% simulate system using these inputs
y0 = 0;
x0 = [0 0]';
ytrain = zeros( numtrials , kfinal );
for i = 1 : numtrials
    [ yout , xout ] = sim_discrete( @vf_siso , kfinal , x0 , y0 , utrain(i,:) );
    ytrain(i,:) = yout';
%     xtrain(i,:) = xout';
end
% ytrain = ytrain(:,2:end);   % remove initial point, it's zero by construction DEBUG
% utrain = utrain(:,1:end-1); % remove last point, it is never actually used DEBUG
% ktrain = 1 : kfinal-1;  % since we removed a pont we have to make this smaller DEBUG

%% simulate bilinear system

xbl0 = zeros( size(A,1) , 1 );
ybl0 = 0;

ybl = zeros( kfinal , length(ybl0) );
xbl = zeros( kfinal , length(xbl0) );

ybl(1,:) = ybl0';
xbl(1,:) = xbl0';
for i = 2 : kfinal
    [ yplus , xplus ] = vf_bilinear( (i-1) , xbl(i-1,:)' , utrain , A , N , B , C );
    ybl(i,:) = yplus';
    xbl(i,:) = xplus';
end


end