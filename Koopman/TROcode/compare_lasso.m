% compare_lasso.m
%
% Compare the performance of models across lasso parameters

%% load candidate models

% Class object with no models trained
temp = load('systems\waves\old-robot_class-nomodels.mat');
class_nomodels = temp.sysid_class;

% LASSO with projection
temp = load('systems\waves\old-robot_lasso-candidate-models.mat');
models_proj = temp.candidates;

ksysid_proj = class_nomodels;
ksysid_proj.candidates = models_proj;

% LASSO without projectio
temp = load('systems\waves\old-robot_lasso-noproj-candidate-models.mat');
models_noproj = temp.candidates;

ksysid_noproj = class_nomodels;
ksysid_noproj.candidates = models_noproj;

% LASSO wieght parameters
lasso_params = 0.05 : 0.5 :10.05;    % lasso weights used to train models

%% Run validation on all the models

% LASSO with projection:
% run validation trials for the lasso+projection model
res_proj = cell( size(ksysid_proj.candidates) );    % store results in a cell array
err_proj = cell( size(ksysid_proj.candidates) );    % store error in a cell array 
if iscell(ksysid_proj.candidates)
    for i = 1 : length(ksysid_proj.candidates)
        [ res_proj{i} , err_proj{i} ] = ksysid_proj.valNplot_model( i , false );
    end
else
    [ res_proj{1} , err_proj{1} ] = ksysid_proj.valNplot_model;
end

% compile the average model error as a function of lasso parameter
err_proj_v_lasso = zeros( size(lasso_params) );
for i = 1 : length( err_proj )
    for j = 1 : length( err_proj{i} )
        err_proj_v_lasso(i) = err_proj_v_lasso(i) + err_proj{i}{j}.rmse2;
    end
    err_proj_v_lasso(i) = err_proj_v_lasso(i) / j;  % take average over val trials
end


% LASSO without projection:
% run validation trials for the lasso+projection model
res_noproj = cell( size(ksysid_noproj.candidates) );    % store results in a cell array
err_noproj = cell( size(ksysid_noproj.candidates) );    % store error in a cell array 
if iscell(ksysid_noproj.candidates)
    for i = 1 : length(ksysid_noproj.candidates)
        [ res_noproj{i} , err_noproj{i} ] = ksysid_noproj.valNplot_model( i , false );
    end
else
    [ res_noproj{1} , err_noproj{1} ] = ksysid_noproj.valNplot_model;
end

% compile the average model error as a function of lasso parameter
err_noproj_v_lasso = zeros( size(lasso_params) );
for i = 1 : length( err_noproj )
    for j = 1 : length( err_noproj{i} )
        err_noproj_v_lasso(i) = err_noproj_v_lasso(i) + err_noproj{i}{j}.rmse2;
    end
    err_noproj_v_lasso(i) = err_noproj_v_lasso(i) / j;  % take average over val trials
end

%% smooth out the erro data a little bit

% remove outliers
err_proj_v_lasso_noout = filloutliers(err_proj_v_lasso,'linear');
err_noproj_v_lasso_noout = filloutliers(err_noproj_v_lasso,'linear');

% take moving average
err_proj_v_lasso_smooth = movmean( err_proj_v_lasso_noout , 8 );
err_noproj_v_lasso_smooth = movmean( err_noproj_v_lasso_noout , 8 );


%% Determine the A matrices of the models

dens_proj = zeros( size(lasso_params) );
dens_noproj = zeros( size(lasso_params) );
for i = 1 : length( lasso_params )
    dens_proj(i) = nnz( ksysid_proj.candidates{i}.A( find( abs( ksysid_proj.candidates{i}.A ) > 5e-6 ) ) ) / numel( ksysid_proj.candidates{i}.A );
    dens_noproj(i) = nnz( ksysid_noproj.candidates{i}.A( find( abs( ksysid_noproj.candidates{i}.A ) > 5e-6 ) ) ) / numel( ksysid_noproj.candidates{i}.A );
end

%% plot it

figure; 
hold on;
yyaxis left;
plot( lasso_params(1:10)  , dens_proj(10:-1:1) , 'r');
plot( lasso_params(1:10)  , dens_noproj(10:-1:1) , 'b' );
yyaxis right;
semilogy( lasso_params(1:10) , err_proj_v_lasso_smooth(1:10) , '-r');
semilogy( lasso_params(1:10) , err_noproj_v_lasso_smooth(1:10) , '-b');








        