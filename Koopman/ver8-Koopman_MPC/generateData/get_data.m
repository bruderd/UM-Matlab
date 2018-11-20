function [ trialData ] = get_data( data_file, data_path, params )
%get_data: Read in empirical data with time, state, and input information
%   -Valid input file type is .mat
%   -Assumes you load in a single time series
%   -It interpolates input data so that data points are spaced by
%   sample period Ts.

numericalDerivs = params.numericalDerivs;

% create data struct
trialData = struct;

[data_filepath,data_name,data_ext] = fileparts([data_path, data_file]);

%% Read in data file, and interpolate so points are Ts apart
if data_ext == '.mat'
    
    raw = load([data_path, data_file]);
    traw = raw.t;
    
%   % STRETCH OUT TIME
%     traw = 12 * traw;    
    
    traw = traw - traw(1);  % remove any offset in start time
    xraw = raw.x;
    uraw = raw.u;
    [traw_uq, index] = unique(traw);    % filter out non-unique values
    
    if numericalDerivs    % do this if numerical derivatives are req.
        % filter raw data to smooth it out a bit
        xraw = movmean(xraw, 50);
        
        % interpolate data so that it lines up with sampling times
        tq = ( 0:params.Ts:traw(end) )';
        xq = interp1(traw_uq, xraw(index,:), tq, 'spline');   % interpolate results to get samples at sampling interval Ts
        uq = interp1(traw_uq, uraw(index,:), tq, 'spline');
        
        % filter state measurements to lessen noise impact
        xfilt = movmean(xq, params.filterWindow(1));
        
%         % ISOLATE xy coordinates
%         xfilt = xfilt(:,1:2);
        
        % take numerical derivatives between sampled points
        numstates = params.n/2;
        xdot = ( xfilt(2:end, :) - xfilt(1:end-1,:) ) / params.Ts;
        xdot = [xdot; zeros(1,numstates)];  % pad end with zeros to make size consistent
        
        % filter numerical derivatives to lessen noise impact
        xdotfilt = movmean(xdot, params.filterWindow(2));
        x = [xfilt, xdotfilt];
        
%         % SCALE DATA (REMOVE THIS LATER!!)
%         [x, uq] = scale_data(x, uq, params);
        
        % define output
        trialData.t = tq;
        trialData.x = x;
        trialData.u = uq;
        trialData.y = trialData.x;    % observed state == state since we have no other ground truth
    else 
        % filter raw data to smooth it out a bit
        xraw = movmean(xraw, 50);     
        
        % interpolate data so that it lines of with sampling times
        tq = ( 0:params.Ts:traw(end) )';
        xq = interp1(traw,xraw,tq, 'spline');   % interpolate results to get samples at sampling interval Ts
        uq = interp1(traw,uraw,tq, 'spline');
        
        % filter state measurements to lessen noise impact
        xfilt = movmean(xq, params.filterWindow(1));
        
%         % SCALE DATA (REMOVE THIS LATER!!)
%         [xfilt, uq] = scale_data(xfilt, uq, params);

%         % ISOLATE xy coordinates
%         xfilt = xfilt(:,1:2);

        % define output
        trialData.t = tq;
        trialData.x = xfilt;
        trialData.u = uq;
        trialData.y = trialData.x;    % observed state == state since we have no other ground truth
    end
    
else 
    
    disp('Invalid file type selected. Data must be in .mat format')
    
end

end

