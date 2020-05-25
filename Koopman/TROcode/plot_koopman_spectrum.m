% plot_koopman_spectrum
%
% plot the eigenvalues of the discrete model A matrix found using the
% Koopman operator on the complex plane


%% get eigenvalues of the A matrix

evals = eig( A );
evals_real = real( evals );
evals_imag = imag( evals );

%% get points to plot the unit circle

t = 0 : 0.01 : 2*pi + 0.01;
circ_real = cos( t );
circ_im = sin( t );

figure;
xlim([-1.5 , 1.5]);
ylim([-1.5,1.5]);
xL = xlim;
yL = ylim;
hold on;
line( [0 0] , yL , 'Color' , 'k');  %y-axis
line( xL , [0 0] , 'Color' , 'k');  %x-axis
plot( circ_real , circ_im , 'k--');    % unit circle
plot( evals_real , evals_imag , '*');  % eigenvalues
hold off;
grid on;
daspect([1 1 1]);   % make axis ratio 1:1

% zoomed in
figure;
xlim([0.92 1.02]);
ylim([-0.05 0.05]);
xL = xlim;
yL = ylim;
hold on;
line( xL , [0 0] , 'Color' , 'k');  %x-axis
plot( circ_real , circ_im , 'k--');    % unit circle
plot( evals_real , evals_imag , '*');  % eigenvalues
hold off;
box on;
% daspect([1 1 1]);   % make axis ratio 1:1