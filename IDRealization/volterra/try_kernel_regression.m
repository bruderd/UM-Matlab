% try_kernel_regression

%% simulate system

kfinal = 1000;

% generate input as a random walk in [-1,1]
steps = randi([-1 1],1,kfinal);
scales = rand(1,kfinal);
u = zeros( [kfinal , 1] );
for i = 2 : length(steps)
    if u(i-1) >= 1
        u(i) = u(i-1) - 1e-1;
    elseif u(i-1) <= -1
        u(i) = u(i-1) + 1e-1;
    else
        u(i) = u(i-1) + 2e-1 * scales(i) * steps(i);
    end
end

% simulate system using these inputs
y0 = 0;
x0 = [0 0]';
[ yout , xout ] = sim_discrete( @vf_siso , kfinal , x0 , y0 , u );
uout = u;
kout = 1 : kfinal;

%% identify volterra series

hor = 10; % length of model horizon
p = 1; % order of discrete volterra series

% partition data
K = hor : hor : kfinal;
N = length(K); % number of training points (wrt the output)
U = zeros( length(K) , hor );
for i = 1 : length(K)
    U(i,:) = uout( (i-1)*hor+1 : i*hor )';
end
Y = yout(K);

% build data matrix for regression
data_output = Y;
data_input = zeros(N,N);
for i = 1 : N
    for j = 1 : N
        data_input(i,j) = ( 1 + U(j,:) * U(i,:)' )^p;
    end
end

% solve for the coefficients
Alpha = lsqminnorm( data_input , data_output );

%% validate the identified model

% generate new set of inputs as a random walk in [-1,1]
uval = zeros( [kfinal , 1] );
for i = 2 : length(steps)
    if uval(i-1) >= 1
        uval(i) = uval(i-1) - 1e-1;
    elseif uval(i-1) <= -1
        uval(i) = uval(i-1) + 1e-1;
    else
        uval(i) = uval(i-1) + 2e-1 * scales(i) * steps(i);
    end
end

% simulate using the real model
[ yout_real , xout_real ] = sim_discrete( @vf_siso , kfinal , x0 , y0 , uval );
uout_real = uval;

% simulate using the identified model
yout_fake = zeros( kfinal , 1 );
uout_fake = [ zeros(hor,1) ; uout_real ];
for i = 2 : kfinal
    yout_fake(i) = Alpha' * ( 1 + U * uout_fake(i:i+hor-1) ).^p;
end

% plot the real model verse the fake one
figure;
hold on;
plot(yout_real);
plot(yout_fake);
hold off;
legend('Real','Volterra');









