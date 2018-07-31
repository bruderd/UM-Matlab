function [ data ] = get_data( params, numericalDerivs )
%get_data: Read in empirical data with time, state, and input information
%   -Valid input file types are .csv and .mat
%   -Assumes you load in a single time series
%   -It interpolates input data so that data points are spaced by
%   sample period Ts.

% check if the user state includes derivatives or if numerical derivatives
% are required
if (~exist('numericalDerivs', 'var'))
    numericalDerivs = 'off';
end

% create data struct
data = struct;

% Prompt user to identify data file
[data_file,data_path] = uigetfile;
[data_filepath,data_name,data_ext] = fileparts([data_path, data_file]);

% Read in data file
if data_ext == '.csv'
    
    if numericalDerivs == 'on'  % do this if numerical derivatives are req.
        M = csvread([data_path, data_file]);
        traw = M(:, 1);
        xraw = M(:, 2 : 2+(params.n/2));
        uraw = M(:, 2+(params.n/2) : 2+(params.n/2)+params.p);
        
        % interpolate data so that it lines of with sampling times
        tq = 0:params.Ts:traw(end);
        xq = interp1(traw,xraw,tq);   % interpolate results to get samples at sampling interval Ts
        uq = interp1(traw,uraw,tq);
        
        % calculate numerical derivative at each point
        xqdot = ( xq(2:end,:) - xq(1:end-1) ) / params.Ts;
        xqdot = [xqdot; zeros(1,params.n/2)];
        
        % define output
        data.t = tq;
        data.x = [xq, xqdot];
        data.u = uq;
    else
        M = csvread([data_path, data_file]);
        data.t = M(:, 1);
        data.x = M(:, 2 : 2+params.n);
        data.u = M(:, 2+params.n : 2+params.n+params.p);
        
        % interpolate data so that it lines of with sampling times
        tq = 0:params.Ts:traw(end);
        xq = interp1(traw,xraw,tq);   % interpolate results to get samples at sampling interval Ts
        uq = interp1(traw,uraw,tq);

        % define output
        data.t = tq;
        data.x = xq;
        data.u = uq;
    end
    
elseif data_ext == '.mat'
    
    raw = load([data_path, data_file]);
    traw = raw.t;
    xraw = raw.x;
    uraw = raw.u;
    if numericalDerivs == 'on'  % do this if numerical derivatives are req.
        % interpolate data so that it lines of with sampling times
        tq = 0:params.Ts:traw(end);
        xq = interp1(traw,xraw,tq);   % interpolate results to get samples at sampling interval Ts
        uq = interp1(traw,uraw,tq);
        
        % calculate numerical derivative at each point
        xqdot = ( xq(2:end,:) - xq(1:end-1) ) / params.Ts;
        xqdot = [xqdot; zeros(1,params.n/2)];
        
        % define output
        data.t = tq;
        data.x = [xq, xqdot];
        data.u = uq;
    else 
        % interpolate data so that it lines of with sampling times
        tq = 0:params.Ts:traw(end);
        xq = interp1(traw,xraw,tq);   % interpolate results to get samples at sampling interval Ts
        uq = interp1(traw,uraw,tq);

        % define output
        data.t = tq;
        data.x = xq;
        data.u = uq;
    end
    
else 
    
    disp('Invalid file type selected. Data must be in .csv or .mat format')
    
end

end
