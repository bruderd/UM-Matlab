function [ x_scaled, u_scaled , w_scaled , params ] = scale_data( x , u , w , params)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

scale = params.scale;
n = size(x,2);
p = size(u,2);
nw = size(w,2);

% maximum value for each of the states (doesn't include input), in a row vector
maxStates = max( abs( max(x, [],1) ) , abs( min(x, [],1) ));
maxInput = max( abs( max(u, [],1) ) , abs( min(u, [],1) ));
maxLoad = max( abs( max(w, [],1) ) , abs( min(w, [],1) ));

% scale each state so that it never goes outside the range [-scale, scale]
x_scaled = scale * x ./ maxStates;
x_scaled(isnan(x_scaled))=0;    % convert any NaNs to zeros

% scale the input so that it never goes above scale
u_scaled = scale * u ./ maxInput;
u_scaled(isnan(u_scaled))=0;    % convert any NaNs to zeros

% scale the input so that it never goes above scale
w_scaled = scale * w ./ maxInput;
w_scaled(isnan(w_scaled))=0;    % convert any NaNs to zeros

%% store the scaling factor in the params struct
params.xScaleFactor = scale ./ maxStates;
params.uScaleFactor = scale ./ maxInput;
params.wScaleFactor = scale ./ maxLoad;

% convert any infs to zero
params.xScaleFactor(isinf(params.xScaleFactor))=1;    % convert any Infs to ones
params.uScaleFactor(isinf(params.uScaleFactor))=1;    % convert any Infs to ones
params.wScaleFactor(isinf(params.wScaleFactor))=1;    % convert any Infs to ones

end

