function plot_state_eval( data , koopsim , params)
%plot_state_eval: Creates plot that shows the agreement of the
%state predictions given a specific Koopman operator
%   Detailed explanation goes here

fig = figure;

range = 1 : 12*150;  % data index range 

% scale the data to be in centimeters (TO BE COMPLETED...)
% valdata_cm = data.val1.x * ( params.

% f(x) = x_1
subplot(2,1,1)
hold on
plot( data.val1.t(range) , data.val1.x(range,1) , 'Color', [0.5 , 0.5 , 0.5] , 'LineWidth' , 2);  % real data
plot( koopsim.val1.t(range) , koopsim.val1.x(range,1) , 'b:' , 'LineWidth' , 2);  % koopman simulation
hold off
ylabel( '$ y $ (cm)' , 'interpreter' , 'latex' );
xlabel( 'Time (s)' , 'interpreter' , 'latex' );
yticks([]);

% create legend for the plots
legend({ 'Validation Data' , 'Koopman Prediction' } , 'Location' , 'northeast');

% f(x) = x_2
subplot(2,1,2)
hold on
plot( data.val1.t(range) , data.val1.x(range,2) , 'Color', [0.5 , 0.5 , 0.5] , 'LineWidth' , 2);  % real data
plot( koopsim.val1.t(range) , koopsim.val1.x(range,2) , 'b:' , 'LineWidth' , 2);  % koopman simulation
hold off
ylabel( '$ z $ (cm)' , 'interpreter' , 'latex' );
xlabel( 'Time (s)' , 'interpreter' , 'latex' );
yticks([]);


% % change y-axis of all plots in figure
% fig = gcf;
% allaxes = findall(fig, 'type', 'axes');
% allaxes(5).YLim = [-params.scale , params.scale];
% allaxes(6).YLim = [-params.scale , params.scale];
        
end