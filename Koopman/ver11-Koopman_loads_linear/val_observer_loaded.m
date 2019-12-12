function [error, koopsim] = val_observer_loaded( data, lifted )
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
liftW = str2func([ 'Wlift_' , params.systemName ]);
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
    psi0 = liftState( zeta0 , valdata.w(1,:)' );
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
            psireliftk = liftState(zrelift , valdata.w(i,:)');
        end
        % linear dynamics
        psikp1 = lifted.A * psik + lifted.B * valdata.u(i,:)';
        xdis(i+1,:) = ( lifted.C * psikp1 )';
        
        % dynamics when relifting at each timestep
        psireliftkp1 = lifted.A * psireliftk + lifted.B * valdata.u(i,:)';
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
    what = zeros( length(tspan) , params.nw);
    whatkp1 = zeros( params.nw , 1 );
    
    % initialize the state with load appended at the end
    xhat = [ y , what ];
    
    % build matrices for estimating the load
    hor = 20;    % number of time steps to be considered in load estimate equation
    UK = zeros(params.p*hor,1); % stack of inputs
    Ykp1 = zeros(params.ny*hor,1);
    Wk = zeros(params.N*hor , 1 + params.nw);
    Astack = kron(eye(hor) , lifted.A);
    Bstack = kron(eye(hor) , lifted.B);
    CAstack = kron(eye(hor) , lifted.Cy * lifted.A);
    CBstack = kron(eye(hor) , lifted.Cy * lifted.B);
    
    for i = 1 : length(tspan) - 1
        
        % Lift the previous measurement
        if params.nd == 1   % assuming 1 or zeros delays, not 2 or more (for now)
            if i == 1
                zetak = [ y(i,:)' ; y(i,:)' ; valdata.u(i,:)' ]; % the nonlifted state when there is 1 delay
            else
                zetak = [ y(i,:)' ; y(i-1,:)' ; valdata.u(i,:)' ]; % the nonlifted state when there is 1 delay
            end
        else
            zetak = y(i,:)';
        end
        Whatk = liftW(zetak);
        psihat(i,:) = Whatk * [1 ; whatkp1]; % used previous estimate of load to update psihat
        
        % Measurement
        mu = 0; % mean of measurement noise
        sigma = 0.00;  % standard deviation of measurement noise
%         ykp1 = xdis(i+1,1:params.ny)' + normrnd(mu,sigma);   % has some noise added in
        ykp1 = xreal(i+1,1:params.ny)' + normrnd(mu,sigma);   % measure the real nonl
%         ykp1 = psirelift(i+1,1:params.ny)' + normrnd(mu,sigma);   % has some noise added in
        y(i+1,:) = ykp1';
        
%         % Lift the measured output (should do this before measuring)
%         if params.nd == 1   % assuming 1 or zeros delays, not 2 or more (for now)
%             zetakp1 = [ ykp1 ; y(i,:)' ; valdata.u(i,:)' ]; % the nonlifted state when there is 1 delay
%         else
%             zetakp1 = ykp1;
%         end
%         Whatkp1 = liftW(zetakp1);
%         psihat(i+1,:) = Whatkp1 * [1 ; whatkp1]; % used previous estimate of load to update psihat
        
        % get estimate of the load using least squares over past horizon (on xhat)
        Ykp1 = [ ykp1 ; Ykp1(1 : params.ny*(hor-1)) ];
        Wk = [ Whatk ; Wk(1 : params.N*(hor - 1) , :) ];
        UK = [ valdata.u(i,:)' ; UK(1 : params.p*(hor-1)) ];
        Clsqlin = CAstack * Wk;
        dlsqlin = Ykp1 - CBstack * UK;
        Aeq = blkdiag( 1 , zeros(params.nw , params.nw) );
        beq = [ 1 ; zeros(params.nw,1) ];
        lb = zeros(params.nw+1,1);
        ub = [ Inf ; 0.9 * ones(params.nw,1) ];
        sol = lsqlin( Clsqlin , dlsqlin , [] , [] , Aeq , beq , lb , ub );  % solve for what using constrained least squares solver
        whatkp1 = sol(2:end);
%         if all( whatkp1 < ones(params.nw,1) ) && all( whatkp1 >= zeros(params.nw,1) )  % set what to zero if solution lies outside of bounds [0,1]
%             what(i+1,:) = whatkp1';
%         else
%             what(i+1,:) = zeros(1,params.nw);
%         end
        
%         % Use lasso instead of just least squares
%         lasso = 0; % L1 penalty weight
%         H = 2 * ( CWstack' * CWstack + lasso );
%         f = -2 * CWstack' * ( YKp1 - CAstack * PsiK - CBstack * UK );
%         whatkp1 = quadprog(H,f,[],[],[],[],0,0.9);
%         what(i+1,:) = whatkp1';
       

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