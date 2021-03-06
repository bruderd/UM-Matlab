function plot_comparison_all_models( err , res , ind )
%Plots the model predictions found in 'evaluate_all_models'
%   Does the same thing as plot_modelPrediction_wNLmodel_wNnet but with
%   data other than the sinusoidal noise data.
%   
%   err: output from 'evaluate_all_models.m'
%   res: output from 'evaluate_all_models.m' 
%   ind: index of the experiment you want to generate the plot from
%   
%   FIXED: This did not include a conversion from inches to cm until
%   2019-11-11. Now it is included but only for layout version 3.

% Define colors
cb_red = [228,26,28] ./ 255;   % red
cb_blue = [55,126,184] ./ 255;  % blue
cb_orange = [255,127,0] ./ 255;   % orange
cb_brown = [166,86,40] ./ 255;  % brown


%% LAYOUT FROM PAPER
fig1 = figure;

% subplot(3,2,1)
% hold on;
% plot( in2cm( res.real{ind}.y(:,1) ) , in2cm( res.real{ind}.y(:,2) ) , 'Color' , [0 0 0 1] ,  'LineWidth' , 2 , 'LineStyle' , '-' );
% plot( in2cm( res.kooplin{ind}.y(:,1) ) , in2cm( res.kooplin{ind}.y(:,2) ) , 'Color' , cb_blue ,  'LineWidth' , 2 , 'LineStyle' , '--');
% plot( in2cm( data{ind}.horizon.Ylin(:,1) ) , in2cm( data{ind}.horizon.Ylin(:,2) ) , 'Color' , cb_red ,  'LineWidth' , 2 , 'LineStyle' , ':');
% plot( in2cm( data{ind}.horizon.Ylin(:,1) ) , in2cm( data{ind}.horizon.YNLkoop(:,2) ) , 'Color' , cb_orange ,  'LineWidth' , 2 , 'LineStyle' , '-.');
% plot( in2cm( data{ind}.horizon.Ynnet(:,1) ) , in2cm( data{ind}.horizon.Ynnet(:,2) ) , 'Color' , cb_brown ,  'LineWidth' , 2 , 'LineStyle' , ':');
% hold off;
% ylabel('$x_2$ (cm)' , 'Interpreter' , 'Latex' , 'FontSize', 16)
% xlabel('$x_1$ (cm)' , 'Interpreter' , 'Latex' , 'FontSize', 16)
% % xlabel('Time (s)')
% ylim( in2cm( data{ind}.horizon.Yreal(1,2) ) + [-15 , 15] );
% % xlim([0 , data{ind}.horizon.T(end)]);
% xlim( in2cm( data{ind}.horizon.Yreal(1,1) ) + [-15 , 15] );
% grid on;
% box on;

x_limits = [ -10 , 15 ];
y_limits = [ -10 , 15 ];

% Linear Koopman model
subplot(3,2,1)
title('Koopman (linear)', 'FontSize', 12);
hold on;
plot( in2cm( res.real{ind}.y(:,1) ) , in2cm( res.real{ind}.y(:,2) ) , 'Color' , [0.5 0.5 0.5 1] ,  'LineWidth' , 2 , 'LineStyle' , '-' );
plot( in2cm( res.real{ind}.y(end,1) ) , in2cm( res.real{ind}.y(end,2) ) , '.' , 'Color' , [0.5 0.5 0.5 1] , 'MarkerSize' , 16);
plot( in2cm( res.kooplin{ind}.y(:,1) ) , in2cm( res.kooplin{ind}.y(:,2) ) , 'Color' , cb_blue ,  'LineWidth' , 2 , 'LineStyle' , '--');
plot( in2cm( res.kooplin{ind}.y(end,1) ) , in2cm( res.kooplin{ind}.y(end,2) ) , '.' , 'Color' , cb_blue , 'MarkerSize' , 16);
hold off;
ylabel('$x_2$ (cm)' , 'Interpreter' , 'Latex' , 'FontSize', 16)
% xlabel('$x_1$ (cm)' , 'Interpreter' , 'Latex' , 'FontSize', 16)
ylim( y_limits );
xlim( x_limits );
grid on;
box on;

% Linear model
subplot(3,2,2)
title('State Space', 'FontSize', 12);
hold on;
plot( in2cm( res.real{ind}.y(:,1) ) , in2cm( res.real{ind}.y(:,2) ) , 'Color' , [0.5 0.5 0.5 1] ,  'LineWidth' , 2 , 'LineStyle' , '-' );
plot( in2cm( res.real{ind}.y(end,1) ) , in2cm( res.real{ind}.y(end,2) ) , '.' , 'Color' , [0.5 0.5 0.5 1] , 'MarkerSize' , 16);
plot( in2cm( res.lin{ind}.y(:,1) ) , in2cm( res.lin{ind}.y(:,2) ) , 'Color' , cb_red ,  'LineWidth' , 2 , 'LineStyle' , ':');
plot( in2cm( res.lin{ind}.y(end,1) ) , in2cm( res.lin{ind}.y(end,2) ) , '.' , 'Color' , cb_red , 'MarkerSize' , 16);
hold off;
% ylabel('$x_2$ (cm)' , 'Interpreter' , 'Latex' , 'FontSize', 16)
% xlabel('$x_1$ (cm)' , 'Interpreter' , 'Latex' , 'FontSize', 16)
ylim( y_limits );
xlim( x_limits );
grid on;
box on;

% Nonlinear Koopman model
subplot(3,2,3)
title('Koopman (nonlinear)', 'FontSize', 12);
hold on;
plot( in2cm( res.real{ind}.y(:,1) ) , in2cm( res.real{ind}.y(:,2) ) , 'Color' , [0.5 0.5 0.5 1] ,  'LineWidth' , 2 , 'LineStyle' , '-' );
plot( in2cm( res.real{ind}.y(end,1) ) , in2cm( res.real{ind}.y(end,2) ) , '.' , 'Color' , [0.5 0.5 0.5 1] , 'MarkerSize' , 16);
plot( in2cm( res.koop{ind}.y(:,1) ) , in2cm( res.koop{ind}.y(:,2) ) , 'Color' , cb_orange ,  'LineWidth' , 2 , 'LineStyle' , '-.');
plot( in2cm( res.koop{ind}.y(end,1) ) , in2cm( res.koop{ind}.y(end,2) ) , '.' , 'Color' , cb_orange , 'MarkerSize' , 16);
hold off;
ylabel('$x_2$ (cm)' , 'Interpreter' , 'Latex' , 'FontSize', 16)
xlabel('$x_1$ (cm)' , 'Interpreter' , 'Latex' , 'FontSize', 16)
ylim( y_limits );
xlim( x_limits );
grid on;
box on;

% Nonlinear Koopman model
subplot(3,2,4)
title('Neural Network', 'FontSize', 12);
hold on;
plot( in2cm( res.real{ind}.y(:,1) ) , in2cm( res.real{ind}.y(:,2) ) , 'Color' , [0.5 0.5 0.5 1] ,  'LineWidth' , 2 , 'LineStyle' , '-' );
plot( in2cm( res.real{ind}.y(end,1) ) , in2cm( res.real{ind}.y(end,2) ) , '.' , 'Color' , [0.5 0.5 0.5 1] , 'MarkerSize' , 16);
plot( in2cm( res.nnet{ind}.y(:,1) ) , in2cm( res.nnet{ind}.y(:,2) ) , 'Color' , cb_brown ,  'LineWidth' , 2 , 'LineStyle' , ':');
plot( in2cm( res.nnet{ind}.y(end,1) ) , in2cm( res.nnet{ind}.y(end,2) ) , '.' , 'Color' , cb_brown , 'MarkerSize' , 16);
hold off;
% ylabel('$x_2$ (cm)' , 'Interpreter' , 'Latex' , 'FontSize', 16)
xlabel('$x_1$ (cm)' , 'Interpreter' , 'Latex' , 'FontSize', 16)
ylim( y_limits );
xlim( x_limits );
grid on;
box on;


% error plot
subplot(3,2,[5,6])
hold on;
p1 = plot( res.real{ind}.t , zeros( size( res.kooplin{ind}.err ) ), 'Color' , [0.5 0.5 0.5 1] ,  'LineWidth' , 2 );
p2 = plot( res.real{ind}.t , in2cm( res.kooplin{ind}.err ) , 'Color' , cb_blue ,  'LineWidth' , 2 , 'LineStyle' , '--');
p3 = plot( res.real{ind}.t , in2cm( res.koop{ind}.err ) , 'Color' , cb_orange ,  'LineWidth' , 2 , 'LineStyle' , '-.');
p4 = plot( res.real{ind}.t , in2cm( res.lin{ind}.err ) , 'Color' , cb_red ,  'LineWidth' , 2 , 'LineStyle' , ':');
p5 = plot( res.real{ind}.t , in2cm( res.nnet{ind}.err ) , 'Color' , cb_brown ,  'LineWidth' , 2 , 'LineStyle' , ':');
hold off;
ylabel('Error (cm)' , 'FontSize', 14)
xlabel('Time (s)' , 'FontSize', 14)
legend([p1 , p4 , p2 , p3 , p5] , {'Actual' , 'Linear State-space' , 'Linear Koopman' , 'Nonlinear Koopman' , 'Neural Network'} , 'Location' , 'northeast', 'FontSize', 11);
ylim( [0 , 15] );
xlim( [0 , res.real{ind}.t(end)] );
grid on;
box on;