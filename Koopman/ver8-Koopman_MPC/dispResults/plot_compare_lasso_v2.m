% plot_compare_lasso_v2
%   Generates the figure that shows how sparsity and error are affected by 
%   the choice of lasso parameter

% load the datafile (update this if new comparison file is desired)
load(['lassoCompareData',filesep,'realdata_poly4.mat']);

%% Mess with the data to get it in correct plotting form

% % reverse the direction of the data
% derror = flipud(totnerror.dist);
% Adensity = fliplr(density.A);
% MAdensity = fliplr(density.MA);
% lambda = fliplr(lassoParams);

% assign new names (only needed since above is commented out)
derror = totnerror.dist;
Adensity = density.A;
MAdensity = density.MA;
lambda = lassoParams;

% take every other data point
derror = [ derror(1) ; derror(2:2:end) ];
Adensity = [ Adensity(1) , Adensity(2:2:end) ];
MAdensity = [ MAdensity(1) , MAdensity(2:2:end) ];
lambda = [ 0.01 , lambda(2:2:end) ];    % replace first 0 with 0.01 so log plot can handle it

%% generate figure

% plot the stuff
figure;
semilogx( lambda , derror , 'LineWidth' , 3 , 'LineStyle' , '-' , 'Color' , [117,112,179]/255);
hold on;
% semilogx( lambda , Adensity , 'LineWidth' , 2 );
semilogx( lambda , MAdensity , 'LineWidth' , 3 , 'LineStyle' , ':' , 'Color' , [217,95,2]/255 );
hold off;
set ( gca, 'xdir', 'reverse' ); % reverse xaxis direction

% plot the best choise of lambda
bestlam = lambda(4); % occurs at the 4th spot
hold on;
line( [bestlam bestlam] , get(gca,'YLim') , 'Color' , [0 0 0 0.25] , 'LineWidth' , 3);
line( [bestlam 5] , [MAdensity(4), MAdensity(4)] , 'Color' , [0 0 0 0.5] , 'LineWidth' , 1 , 'LineStyle' , '--');

% modify figure appearance
axis([ 0.01 , 5 , 0 , 1 ]);
legend( { 'Model Prediction Error' , 'Density of $\hat{A}$' } , 'location' , 'southwest' , 'Interpreter' , 'Latex');
xlabel('$\lambda$' , 'Interpreter' , 'Latex');

set(gca,'TickLabelInterpreter', 'Latex');
xticks([0.01 , 5]);
xticklabels({ 50 , 0 });
% xticklabels({ '$0$ solution' , '$L^2$ solution' });
yticks([0, MAdensity(4) , 1]);
yticklabels([0, round( MAdensity(4) ,2) ,1]);
% xtickangle(50);



% figure;
% hold on;
% plot( lambda , derror );
% plot( lambda , Adensity );
% plot( lambda , MAdensity );
% hold off;