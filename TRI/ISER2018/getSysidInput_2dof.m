function Psysid_V = getSysidInput_2dof( params, testname )
%getTestInput_2dof - creates the control input files for sysid on the 2dof system
%   testname is the name of the folder the csv is to be saved in

Psysid_Pa = calcPressureSysidPoints(params);    % sysid control pressure [Pa]
Psysid_V = Pa2V(Psysid_Pa, params);     % sysid control pressure [V]
Psysid_V = [zeros(size(Psysid_V(:,1))), Psysid_V];  % add column of zeros since only using valves 2-4


%% save these pressures as .csv files

% check for optional argument, if given, save csv files with that name
if exist('testname','var')
    current_folder = cd;
    dir = strcat(current_folder, '\testPoints\', testname);
    mkdir(char(dir));
    
    % save csv file of control inputs
    saveas1 = strcat(dir, '\sysid.csv');
    csvwrite(char(saveas1), Psysid_V, 0, 0);
    
    % create blank csv file of test data. Will be filled in by labview
    saveas2 = strcat(dir, '\sysidData.csv');
    csvwrite(char(saveas2), [], 0, 0)
end

end