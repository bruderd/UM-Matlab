function plot_observable_eval( data , koopsim , params)
%plot_observable_eval: Creates plot that shows the agreement of the
%observables given a specific Koopman operator
%   Detailed explanation goes here

fig = figure;

% f(x) = x_1
subplot(3,1,1)
hold on
plot( data.val1.t , data.val1.x(:,1) , 'Color', [0.6 , 0.6 , 0.6] , 'LineWidth' , 4);  % real data
plot( koopsim.val1.t , koopsim.val1.x(:,1) , 'b:' , 'LineWidth' , 2);  % koopman simulation
hold off
ylabel( '$ f_1 (x) = x_1 $' , 'interpreter' , 'latex' );
xlabel( 'Time (s)' , 'interpreter' , 'latex' );
ylim([-1, 1]);

% create legend for the plots
legend({ 'Validation Data' , 'Koopman Prediction' } , 'Location' , 'northeast');

% f(x) = x_2
subplot(3,1,2)
hold on
plot( data.val1.t , data.val1.x(:,2) , 'Color', [0.6 , 0.6 , 0.6] , 'LineWidth' , 4);  % real data
plot( koopsim.val1.t , koopsim.val1.x(:,2) , 'b:' , 'LineWidth' , 2);  % koopman simulation
hold off
ylabel( '$ f_2 (x) = x_2 $' , 'interpreter' , 'latex' );
xlabel( 'Time (s)' , 'interpreter' , 'latex' );
ylim([-1, 1]);

% f(x) = x_1 x_2
x1x2_real = data.val1.x(:,1) .* data.val1.x(:,2);
x1x2_koop = koopsim.val1.x(:,1) .* koopsim.val1.x(:,2);
subplot(3,1,3)
hold on
plot( data.val1.t , x1x2_real , 'Color', [0.6 , 0.6 , 0.6] , 'LineWidth' , 4);  % real data
plot( koopsim.val1.t , x1x2_koop , 'b:' , 'LineWidth' , 2);  % koopman simulation
hold off
ylabel( '$ f_3 (x) = x_1 x_2 $' , 'interpreter' , 'latex' );
xlabel( 'Time (s)' , 'interpreter' , 'latex' );
ylim([-1, 1]);
% ax = gca;
% ax.YLim = [-params.scale , params.scale];


% % change y-axis of all plots in figure
% fig = gcf;
% allaxes = findall(fig, 'type', 'axes');
% allaxes(5).YLim = [-params.scale , params.scale];
% allaxes(6).YLim = [-params.scale , params.scale];
        
% ther should be a '...' below the third plot
end

