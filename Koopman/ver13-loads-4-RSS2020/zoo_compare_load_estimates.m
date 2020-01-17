% zoo_compare_load_estimates

num = length(ksysid.valdata);

% add zeta field to all the trials
for i = 1 : num
    [ ksysid.valdata{i} , ~ ] = ksysid.get_zeta( ksysid.valdata{i} );
end

% % add zeta field to all the trials
% [ ksysid.valdata{1} , ~ ] = ksysid.get_zeta( ksysid.valdata{1} );
% [ ksysid.valdata{2} , ~ ] = ksysid.get_zeta( ksysid.valdata{2} );
% [ ksysid.valdata{3} , ~ ] = ksysid.get_zeta( ksysid.valdata{3} );
% [ ksysid.valdata{4} , ~ ] = ksysid.get_zeta( ksysid.valdata{4} );
% [ ksysid.valdata{5} , ~ ] = ksysid.get_zeta( ksysid.valdata{5} );
% [ ksysid.valdata{6} , ~ ] = ksysid.get_zeta( ksysid.valdata{6} );
% [ ksysid.valdata{7} , ~ ] = ksysid.get_zeta( ksysid.valdata{7} );
% [ ksysid.valdata{8} , ~ ] = ksysid.get_zeta( ksysid.valdata{8} );
% [ ksysid.valdata{9} , ~ ] = ksysid.get_zeta( ksysid.valdata{9} );
% [ ksysid.valdata{10} , ~ ] = ksysid.get_zeta( ksysid.valdata{10} );
% [ ksysid.valdata{11} , ~ ] = ksysid.get_zeta( ksysid.valdata{11} );
% [ ksysid.valdata{12} , ~ ] = ksysid.get_zeta( ksysid.valdata{12} );
% [ ksysid.valdata{13} , ~ ] = ksysid.get_zeta( ksysid.valdata{13} );
% [ ksysid.valdata{14} , ~ ] = ksysid.get_zeta( ksysid.valdata{14} );
% [ ksysid.valdata{15} , ~ ] = ksysid.get_zeta( ksysid.valdata{15} );
% [ ksysid.valdata{16} , ~ ] = ksysid.get_zeta( ksysid.valdata{16} );
% [ ksysid.valdata{17} , ~ ] = ksysid.get_zeta( ksysid.valdata{17} );
% [ ksysid.valdata{18} , ~ ] = ksysid.get_zeta( ksysid.valdata{18} );
% [ ksysid.valdata{19} , ~ ] = ksysid.get_zeta( ksysid.valdata{19} );
% [ ksysid.valdata{20} , ~ ] = ksysid.get_zeta( ksysid.valdata{20} );
% [ ksysid.valdata{21} , ~ ] = ksysid.get_zeta( ksysid.valdata{21} );
% [ ksysid.valdata{22} , ~ ] = ksysid.get_zeta( ksysid.valdata{22} );

%% calculate estimates using all trial data
what = zeros(num,1);
wreal = zeros(num,1);
for i = 1 : num
    % solve once over whole trial
    what(i) = ksysid.observer_load( ksysid.valdata{i}.zeta(1:end,:) , ksysid.valdata{i}.u(2:end,:) );
    wreal(i) = ksysid.valdata{i}.w(end);
    
%     % iterative version
%     [ what_ , ~ , ~ , ~ ] = ksysid.val_observer_load_sparse( 100 , 60 , ksysid.valdata{i} );
%     what(i) = what_(end,:);
%     wreal(i) = ksysid.valdata{i}.w(end);
end

%% plot comparison
figure; bar( 1:num , ( wreal + 1 ) * 150 )
hold on; bar( 1:num , ( what + 1 ) * 150 , 'BarWidth', 0.3 );
ylabel('Load (g)');
xlabel('Trial number');
legend('Actual','Estimate');