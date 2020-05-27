% compare_lasso_v2.m
%
% Compare the performance of models across lasso parameters

%% load candidate models

% LASSO with projection
temp = load('systems\fromData\models-4-lasso-diagram_2020-05-21_15-44.mat');   % very large file, this will take some time
ksysid_proj = temp.ksysid;

% LASSO without projection
ksysid_noproj = ksysid_proj;
for i = 1 : length( ksysid_proj.candidates )
   ksysid_noproj.candidates{i}.A = ksysid_proj.candidates{i}.A_noproj;
   ksysid_noproj.candidates{i}.B = ksysid_proj.candidates{i}.B_noproj;
end

% LASSO wieght parameters
lasso_params = ksysid_proj.lasso;    % lasso weights used to train models

%% Run validation on all the models

% LASSO with projection:
% run validation trials for the lasso+projection model
res_proj = cell( size(ksysid_proj.candidates) );    % store results in a cell array
err_proj = cell( size(ksysid_proj.candidates) );    % store error in a cell array 
if iscell(ksysid_proj.candidates)
    for i = 1 : length(ksysid_proj.candidates)
        for j = 1 :  length( ksysid_proj.valdata )
            res_proj{i}{j} = ksysid_proj.val_model( ksysid_proj.candidates{i} , ksysid_proj.valdata{j} );
            err_proj{i}{j} = ksysid_proj.get_error( res_proj{i}{j}.sim , res_proj{i}{j}.real );
        end
    end
else
    res_proj{1} = ksysid_proj.val_model;
    err_proj{1} = ksysid_proj.get_error( res_proj{1}.sim , res_proj{1}.real );
end

% compile the average model error as a function of lasso parameter
err_proj_v_lasso = zeros( size(lasso_params) );
for i = 1 : length( err_proj )
    for j = 1 : length( err_proj{i} )
        err_proj_v_lasso(i) = err_proj_v_lasso(i) + err_proj{i}{j}.neuclid;
    end
    err_proj_v_lasso(i) = err_proj_v_lasso(i) / j;  % take average over val trials
end


% LASSO without projection:
% run validation trials for the lasso+projection model
res_noproj = cell( size(ksysid_noproj.candidates) );    % store results in a cell array
err_noproj = cell( size(ksysid_noproj.candidates) );    % store error in a cell array 
if iscell(ksysid_noproj.candidates)
    for i = 1 : length(ksysid_noproj.candidates)
        for j = 1 :  length( ksysid_proj.valdata )
            res_noproj{i}{j} = ksysid_noproj.val_model( ksysid_noproj.candidates{i} , ksysid_noproj.valdata{j} );
            err_noproj{i}{j} = ksysid_noproj.get_error( res_noproj{i}{j}.sim , res_noproj{i}{j}.real );
        end
    end
else
    res_noproj{1} = ksysid_noproj.val_model;
    err_noproj{1} = ksysid_noproj.get_error( res_noproj{1}.sim , res_noproj{1}.real );
end

% compile the average model error as a function of lasso parameter
err_noproj_v_lasso = zeros( size(lasso_params) );
for i = 1 : length( err_noproj )
    for j = 1 : length( err_noproj{i} )
        err_noproj_v_lasso(i) = err_noproj_v_lasso(i) + err_noproj{i}{j}.neuclid;
    end
    err_noproj_v_lasso(i) = err_noproj_v_lasso(i) / j;  % take average over val trials

end

%% smooth out the error data a little bit

% remove outliers
err_proj_v_lasso_noout = filloutliers(err_proj_v_lasso,'linear');
err_noproj_v_lasso_noout = filloutliers(err_noproj_v_lasso,'linear');

% take moving average
err_proj_v_lasso_smooth = movmean( err_proj_v_lasso_noout , 10 );
err_noproj_v_lasso_smooth = movmean( err_noproj_v_lasso_noout , 10 );


%% Determine the density and rank of the A matrices of the models

dens_proj = zeros( size(lasso_params) );
dens_noproj = zeros( size(lasso_params) );
rank_proj = zeros( size(lasso_params) );
rank_noproj = zeros( size(lasso_params) );
for i = 1 : length( lasso_params )
    dens_proj(i) = nnz( ksysid_proj.candidates{i}.A( find( abs( ksysid_proj.candidates{i}.A ) > 1e-12 ) ) ) / numel( ksysid_proj.candidates{i}.A );
    dens_noproj(i) = nnz( ksysid_noproj.candidates{i}.A( find( abs( ksysid_noproj.candidates{i}.A ) > 1e-12 ) ) ) / numel( ksysid_noproj.candidates{i}.A );
    
    rank_proj(i) = rank( ksysid_proj.candidates{i}.A ); % rank of A matrix
    rank_noproj(i) = rank( ksysid_noproj.candidates{i}.A ); % rank of A matrix
end

%% smooth out the density a little bit

% % take moving average
% dens_proj_smooth = movmean( dens_proj , 5 , 'Endpoints' , 'shrink' );
% dens_noproj_smooth = movmean( dens_noproj , 5 , 'Endpoints' , 'shrink' );

% lowpass filter
dens_proj_smooth = lowpass( dens_proj , 0.25 );
dens_noproj_smooth = lowpass( dens_noproj , 0.25 );

% lowpass filter
rank_proj_smooth = lowpass( rank_proj , 0.1 );
rank_noproj_smooth = lowpass( rank_noproj , 0.1 );

% round the rank since it has to be a whole number
rank_proj_whole = floor( rank_proj_smooth );
rank_noproj_whole = floor( rank_noproj_smooth );


%% plot it

offset = 30;    % starting index offset

fig = figure;

left_color = [ 27,158,119 ] / 255;  % green % [217,95,2] / 255; % purple
right_color = [117,112,179] / 255;   % orange
set(fig,'defaultAxesColorOrder',[left_color; right_color]);

hold on;
yyaxis left;
d_noproj = plot( lasso_params  , dens_noproj_smooth(end:-1:1) );
d_proj = plot( lasso_params  , dens_proj_smooth(end:-1:1) );
% d_noproj = stairs( lasso_params(offset:end)  , rank_noproj_whole(end-offset+1:-1:1) );
% d_proj = stairs( lasso_params(offset:end)  , rank_proj_whole(end-offset+1:-1:1) , '-.' );
ylim([0,1]);
% ylim([0,40]);
ylabel('Koopman Matrix Rank')

yyaxis right;
e_noproj = plot( lasso_params(offset:end) , err_noproj_v_lasso(end-offset+1:-1:1) );
e_proj = plot( lasso_params(offset:end) , err_proj_v_lasso(end-offset+1:-1:1) , '-.' );
ylim([0,1]);
% set(gca, 'YScale', 'log')
ylabel('Model Prediction Error')
hold off;

xticks([ lasso_params(offset) , lasso_params(end) ]);
xticklabels({ '0' , '50' });
xlim([ lasso_params(offset) , lasso_params(end) ]);
xlabel('L^1 Penalty Term, \lambda' )

% gridLegend( [d_noproj , e_noproj , d_proj , e_proj] , 4 , {'poop','pee','butt','hole'} , 'location' , 'north' );






