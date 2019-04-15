function [error, koopsim] = val_observer( data, lifted )
%val_lobserver: Validates the observer for the lifted system, which tries
%to estimate the load by linearly incorporating it into the lifted state.
%   At each time step the states are measured, then the load is estimated
%   based on the linear dynamics.

% params = data.valparams;    % model parameters
params = lifted.params;

error = struct;     % error results comparing real and koopman system
koopsim = struct;   % simulation results for koopman system

% identify the lifting function
cd('liftingFunctions');
liftState = str2func([ 'lift_' , params.systemName ]);
cd('..');

for j = 1 : params.numVals
    %% real system
     
    % isolate the jth validation trial
    valID = ['val', num2str(j)];
    valdata = data.(valID);
    
    index0 = params.nd + 1;  % index of the first state
    tspan = valdata.t(index0 : end);    % start simulation late so delays can be taken into account
    [treal, xreal] = deal(tspan, valdata.x(index0 : end , :));

    %% simulate the behavior of the learned lifted system
    
    % set initial condition
    x0 = valdata.x(index0 , :)';
    xd0 = reshape( flipud( valdata.x(1 : params.nd , :) )' , [params.n * params.nd , 1] );
    ud0 = reshape( flipud( valdata.u(1 : params.nd , :) )' , [params.p * params.nd , 1] );
    zeta0 = [x0; xd0; ud0];
    
    % simulate the behavior of the discrete linear system
    xdis = zeros(length(tspan) , params.n);
    xdis(1,:) = valdata.x(index0 , :);
    psi0 = liftState(zeta0);
    psi = zeros(length(tspan) , params.N);
    psi(1,:) = psi0;
    psirelift = zeros(length(tspan) , params.N);
    psirelift(1,:) = psi0;
    for i = 1 : length(tspan)-1
        if i == 1
            psik = psi0;
            psireliftk = psi0;
        else
            psik = psikp1;
            psireliftk = liftState(zrelift);
        end
        % linear dynamics
        psikp1 = lifted.A * psik + lifted.B * valdata.u(i,:)' + lifted.W * valdata.w(i,:)';
        xdis(i+1,:) = ( lifted.C * psikp1 )';
        
        % dynamics when relifting at each timestep
        psireliftkp1 = lifted.A * psireliftk + lifted.B * valdata.u(i,:)' + lifted.W * valdata.w(i,:)';
        zrelift = lifted.Cz * psireliftkp1;
       
        % save psihat so that I can plot it
        psi(i+1,:) = psikp1;
        psirelift(i+1,:) = psireliftkp1';    % the lifted state if you relift at each timestep
    end
    
    %% simulate observer system based on lifting and solving
    
    % initialize the measured output
    y = zeros( length(tspan) , params.ny );
    y(1,:) = psirelift(1,1:params.ny);
    
    % initialize estimate of lifted state
    psihat = zeros( length(tspan) , params.N );
    psihat(1,:) = psirelift(1,:)';
    
    % initialize load estimate
    what = zeros( length(tspan) , params.nw );
    
    % initialize the state with load appended at the end
    xhat = [ y , what ];
    
    % build matrices for estimating the load
    hor = 100;    % number of time steps to be considered in load estimate equation
    UK = zeros(params.p*hor,1); % stack of inputs
    PsiK = zeros(params.N*hor,1);  % stack of lifted states (no load)
    PsiKp1 = zeros(params.N*hor,1);  % stack of lifted states (no load)
    Wstack = kron( ones(hor,1) , lifted.W );    % matrix of stacked load matrices
    Wstackinv = pinv(Wstack);
    Astack = kron(eye(hor) , lifted.A);
    Bstack = kron(eye(hor) , lifted.B);
    
    for i = 1 : length(tspan) - 1
        
        % Measurement
        mu = 0; % mean of measurement noise
        sigma = 0.01;  % standard deviation of measurement noise
%         ykp1 = xdis(i+1,1:params.ny)' + normrnd(mu,sigma);   % has some noise added in
%         ykp1 = xreal(i+1,1:params.ny)' + normrnd(mu,sigma);   % has some noise added in
        ykp1 = psirelift(i+1,1:params.ny)' + normrnd(mu,sigma);   % has some noise added in
        y(i+1,:) = ykp1';
        
        % Lift the measured output
        psihatkp1 = liftState(ykp1);
        psihat(i+1,:) = psihatkp1';
        
        % get better estimate of the load using least squares over past horizon (on psihat)
        PsiK = PsiKp1;
        PsiKp1 = [ psihatkp1(1:params.N) ; PsiKp1(1 : params.N*(hor-1)) ];
        UK = [ valdata.u(i,:)' ; UK(1 : params.p*(hor-1)) ];
        whatkp1 = Wstackinv * ( PsiKp1 - Astack * PsiK - Bstack * UK);
        what(i+1,:) = whatkp1';

        % update the state with the load appended at the end
        xhatkp1 = [ ykp1 ; whatkp1 ];
        xhat(i+1,:) = xhatkp1';
            
    end 
    
    %% define outputs
%     error.(valID).terror = terror;
%     error.(valID).xerror = xerror;
%     error.(valID).RMSE = RMSE;
    
    koopsim.(valID).t = tspan;
%     koopsim.(valID).x = xsysid;       % should make this include both later...
%     koopsim.(valID).x = xss;
    koopsim.(valID).x = xdis;
    koopsim.(valID).xhat = xhat;    % observer output
    koopsim.(valID).xreal = xreal;  % real system behavior
    koopsim.(valID).u = valdata.u(index0 : end , :);
    
    koopsim.(valID).psi = psi;
    koopsim.(valID).psihat = psihat;
    koopsim.(valID).psirelift = psirelift;
    
    koopsim.(valID).what = what;
    koopsim.(valID).w = valdata.w(index0 : end , :);
   
    %% plot the results
    
%     if params.ploton
%         figure
%         subplot(3,1,1)
%         plot(treal, xreal(:,1:ceil(params.n)))
%         title('Real system')
%         subplot(3,1,2)
%         plot(tsysid, xsysid(:,1:ceil(params.n)))
%         title('Identified system')
%         subplot(3,1,3)
%         hold on
%         plot(terror, xerror(:,1:ceil(params.n)))
%         plot(terror, xerrormax * ones(size(terror)), '--')
%         title('Error')
%         hold off
%     end

end

end
