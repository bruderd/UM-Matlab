

%% Use Koopman operator to perform sysid
waitbar(.5,progress,'Performing Koopman based system identification...');



%% error
waitbar(0.75,progress,'Comparing to validation data set...');

% [error, xkoop] = koopmanSimulation( data.validation, params, koopman ); % only uses koopman transpose, no ODE
if strcmp(basis, 'fourier')
    [error, koopsim] = koopmanValidation_fourier( data, params, koopman );
elseif strcmp(basis, 'poly')
    [error, koopsim] = koopmanValidation( data, params, koopman );
end


%% compare koopman results to those from sysid toolbox
waitbar(0.85,progress,'Preparing data for Matlab SysId toolbox...');

% convert data to a format matlabs sysid toolbox can use
[zsysid_merged, zval_merged, zsysid, zval] = prep_iddata(data);

% save in struct for output
data4sysid = struct;
data4sysid.sysid_merged = zsysid_merged;
data4sysid.val_merged = zval_merged;
data4sysid.val = zval;
data4sysid.sysid = zsysid;
for k = 1:params.numVals
    valID = ['val', num2str(k)];
    zID = ['z', num2str(k)];
    data4sysid.valkoop.(zID) = iddata(koopsim.(valID).x, data.(valID).u, data.valparams.Ts, 'Name', 'Koopman');
end

% show comparison of Koopman system verses ground truth
if params.ploton
    for k = 1: params.numVals
        zID = ['z', num2str(k)];
        figure
        compare(data4sysid.val.(zID), data4sysid.valkoop.(zID));
    end
end

waitbar(1,progress,'Done.');
close(progress);