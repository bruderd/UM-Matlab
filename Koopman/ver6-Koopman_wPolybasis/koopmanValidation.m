function [error, xsysid] = koopmanValidation( valdata, valparams )
%koopmanValidation: Summary of this function goes here
%   Detailed explanation goes here

%% simulate the behavior of learned system

tspan = valdata.t;
x0sim = valparams.x0; % same initial state as validation data initial state
[tsysid, xsysid] = ode45(@(t,x) vf_koopman(x, get_u(t, x, valdata, valparams)), tspan, x0sim);
[treal, xreal] = deal(valdata.t, valdata.x);


%% quantify the error between real behavior and simulated behavior

terror = treal;
xerror = abs( xreal - xsysid );
xerrormax = max(max(xerror(:,1:valparams.n/2)));
% xerrormin = min(min(xerror(:,1:valparams.n/2)));
RMSE = sqrt( sum( (xreal - xsysid).^2 ) / length(terror) );


%% plot the results

if valparams.ploton
    figure
    subplot(3,1,1)
    plot(treal, xreal(:,1:valparams.n/2))
    title('Real system')
    subplot(3,1,2)
    plot(tsysid, xsysid(:,1:valparams.n/2))
    title('Identified system')
    subplot(3,1,3)
    hold on
    plot(terror, xerror(:,1:valparams.n/2))
    plot(terror, xerrormax * ones(size(terror)), '--')
    title('Error')
    hold off
end


% % animate the results
% animate_doublePendulum(sol_real, sol_sysid, valparams);

%% Define outputs

error = struct;
error.terror = terror;
error.xerror = xerror;
error.RMSE = RMSE;


end


function u = get_u(t, x, valdata, valparams)
%get_u: Interpolates to estimate the value of the input at a specific t

u = interp1(valdata.t, valdata.u, t);

end


