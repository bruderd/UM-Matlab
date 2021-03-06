% try_kernel_regression

%% simulate system

kfinal = 100;
numtrials = 100;

utrain = zeros( numtrials , kfinal );
for j = 1 : numtrials
    % generate input as a random walk in [-1,1]
    steps = randi([-1 1],1,kfinal);
    scales = rand(1,kfinal);
    u = zeros( 1 , kfinal);
    for i = 2 : length(steps)
        if u(i-1) >= 1
            u(i) = u(i-1) - 1e-1;
        elseif u(i-1) <= -1
            u(i) = u(i-1) + 1e-1;
        else
            u(i) = u(i-1) + 2e-1 * scales(i) * steps(i);
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
ktrain = 1 : kfinal;

%% identify volterra series

hor = 100; % length of model horizon
p = 2; % order of discrete volterra series

% partition data
% K = hor : hor : kfinal;
% N = length(K); % number of training points (wrt the output)
% U = zeros( numtrials , hor );
% for i = 1 : numtrials
%     U(i,:) = utrain( i, : );
% end
N = numtrials;
U = utrain;
Y = ytrain(:,hor);  % may be hor + 1

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

% generate new set of validation inputs as a random walk in [-1,1]
uval = zeros( numtrials , kfinal );
for j = 1 : numtrials
    % generate input as a random walk in [-1,1]
    steps = randi([-1 1],1,kfinal);
    scales = rand(1,kfinal);
    u = zeros( 1 , kfinal);
    for i = 2 : length(steps)
        if u(i-1) >= 1
            u(i) = u(i-1) - 1e-1;
        elseif u(i-1) <= -1
            u(i) = u(i-1) + 1e-1;
        else
            u(i) = u(i-1) + 2e-1 * scales(i) * steps(i);
        end
    end
    uval(j,:) = u;
end

% simulate system using these inputs and real model
y0 = 0;
x0 = [0 0]';
yval = zeros( numtrials , kfinal );
for i = 1 : numtrials
    [ yout , xout ] = sim_discrete( @vf_siso , kfinal , x0 , y0 , uval(i,:) );
    yval(i,:) = yout';
%     xtrain(i,:) = xout';
end
kval = 1 : kfinal;

% see if regression solution works on new input data
% build data matrix for regression
Yval_real = yval(:,hor);
data_input_val = zeros(N,N);
for i = 1 : N
    for j = 1 : N
        data_input_val(i,j) = ( 1 + uval(j,:) * uval(i,:)' )^p;
    end
end
Yval_fake = data_input_val * Alpha;


% % simulate using the identified model (OLD VERSION)
% yout_fake = zeros( kfinal , 1 );
% uout_fake = [ zeros(hor,1) ; uout_real ];
% % u = [ zeros(hor,1) ; u ];
% for i = 2 : kfinal
%     yout_fake(i) = Alpha' * ( 1 + U * uout_fake(i:i+hor-1) ).^p;
% %     yout_fake(i) = Alpha' * ( 1 + U * u(i:i+hor-1) ).^p;
% %     yout_fake(i) = Alpha' * ( 1 + U * U(i,:)' ).^p;
% end

% % plot the real model verse the fake one (OLD VERSION)
% figure;
% hold on;
% plot(yout_real);
% plot(yout_fake);
% hold off;
% legend('Real','Volterra');



%% Just do regular regression

% matrix of exponents. Each row gives exponents for 1 monomial
exponents = zeros(1,hor);
for i = 1:p
   exponents = [exponents; partitions(i, ones(1,hor))]; 
end
M = size( exponents , 1 );  % number of coefficients needed to be solved

% build data matrix for regression
reg_output = Y;
reg_input = zeros(N,M);
for i = 1 : N
    for j = 1 : M
        monomial = prod( U(i,:).^exponents(j,:) );
        reg_input(i,j) = monomial;
    end
end

% solve for the coefficients
Gamma = lsqminnorm( reg_input , reg_output );


%% validate the regular regression

% simulate using the identified model
yout_reg = zeros( kfinal , 1 );
uout_reg = [ zeros(hor,1) ; uout_real ];
for i = 2 : kfinal
    monomial_vec = prod( kron( ones(M,1) , uout_reg(i:i+hor-1)' ) .^ exponents, 2 );
    yout_reg(i) = Gamma' * monomial_vec;
end

% plot the real model verse the fake one
figure;
hold on;
plot(yout_real);
plot(yout_reg);
hold off;
legend('Real','Volterra');























