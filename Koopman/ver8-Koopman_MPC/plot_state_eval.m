function plot_state_eval( data , koopsim , params)
%plot_state_eval: Creates plot that shows the agreement of the
%state predictions given a specific Koopman operator
%   Detailed explanation goes here

fig = figure;

% f(x) = x_1
subplot(2,1,1)
hold on
plot( data.val1.t , data.val1.x(:,1) , 'Color', [0.5 , 0.5 , 0.5] , 'LineWidth' , 2);  % real data
plot( koopsim.val1.t , koopsim.val1.x(:,1) , 'b' , 'LineWidth' , 1);  % koopman simulation
hold off
ylabel( '$ x_1 $' , 'interpreter' , 'latex' );
xlabel( 'Time (s)' , 'interpreter' , 'latex' );

% create legend for the plots
legend({ 'Validation Data' , 'Koopman Prediction' } , 'Location' , 'northeast');

% f(x) = x_2
subplot(2,1,2)
hold on
plot( data.val1.t , data.val1.x(:,2) , 'Color', [0.5 , 0.5 , 0.5] , 'LineWidth' , 2);  % real data
plot( koopsim.val1.t , koopsim.val1.x(:,2) , 'b' , 'LineWidth' , 1);  % koopman simulation
hold off
ylabel( '$ x_2 $' , 'interpreter' , 'latex' );
xlabel( 'Time (s)' , 'interpreter' , 'latex' );



% % change y-axis of all plots in figure
% fig = gcf;
% allaxes = findall(fig, 'type', 'axes');
% allaxes(5).YLim = [-params.scale , params.scale];
% allaxes(6).YLim = [-params.scale , params.scale];
        
end