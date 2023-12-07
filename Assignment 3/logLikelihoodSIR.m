% log-likelihood function for an SIR epidemic.
function [ll] = logLikelihoodSIR(theta,data,N)

    % Define model parameters as a structure (N is fixed for all MCMC runs)
    para = struct('beta',theta(1),'gamma',theta(2),"iota",theta(3),"delta",theta(4),'N',N); 

    % Define initial conditions as a structure
    ICs = struct('S',para.N-theta(3),'I',theta(3),'R',0);

    % Safety first!
    if min(theta)>0 && ICs.S>0
        % Start log-likelihood at zero
        ll=0;
        % Solve ODEs to recover model trajectory
        [Classes] = ODE_SIR_model(para,ICs,max(data.t));
   
        % loop over observations in data
        for i=1:length(data.t)
            % find the part of the output that coincides with the data
            f=find(Classes.t==data.t(i));

            % add term to cummulative log-likelihood

            %lambda is delta*Classes.I(f)/N
            %this line below doesn't work because of the factorial term,
            %giving -infinity
            %ll=ll+log(poisspdf(data.x(i),theta(4)*Classes.I(f)/N));

            %above gives -Inf
            lambda = theta(4)*Classes.I(f);
            %bigg = log(factorial(data.x(i))); redundant term
            ll = ll + data.x(i)*log(lambda) - lambda;
                       
        end
    else
        ll=-Inf;
    end

    
    