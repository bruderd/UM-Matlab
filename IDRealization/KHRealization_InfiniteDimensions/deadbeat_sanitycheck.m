% dead-beat observer sanity check

% define linear system in observable canonical form
A = [2,1 ; -1.0001,0];
B = 1e-3 * [ -0.0991 ; 0.1991 ];
C = [1 0];
D = 0;

% simulate system for a give set of inputs
k = 1:1000;
u = 0.5 .* rand( length( k ), 1 );

x = zeros( length(k) , 2);
y = zeros( length(k) , 1);
x(1,:) = [ 1 , 0 ];   % initial condition
y(1,:) = C * x(1,:)';
for i = 2 : k(end)
    x(i,:) = ( A * x(i-1,:)' + B * u(i-1) )';
    y(i) = C * x(i,:)';
end

% build dead-beat observer and simulate it
L = 1e-2 * A(:,1);
yhat = zeros(1000,1);
xhat = zeros(1000,2);
for i = 2 : k(end)
    xhat(i,:) = ( A * x(i-1,:)' + B * u(i-1) + L * ( y(i-1) - yhat(i-1) ) )';
    yhat(i) = C * xhat(i,:)';
end

% plot results
figure;
hold on;
plot( k , y );
plot( k ,yhat );
hold off;
legend( 'real' , 'observer' );