%mcmc algorithm function

%takes number of iterations, starting guess theta_0, sigma, n0, and data

%returns

function [output] = mcmc(iters, theta_0, sigma_0, n_0, data)




%acceptance rate monitors
accept = 0;
reject = 0;

%put theta_0 into theta (for readability), current value for MCMC is theta
theta = theta_0;
N = 5500;

%allocate storage for output
stored=zeros(iters,length(theta));
incidence=zeros(iters,1+max(data.t));

%init proposal covariance matrix is sigma_0
sigma = sigma_0;
mu = 0;
epsilon = 10^-6;

%loop
for i=1:iters
    
    %generate proposal for theta and check if it's in range
    proposal = mvnrnd(theta,2.38^2/4*sigma + epsilon*eye(4));
    
    %check, then run
    if(min(proposal)>0 && proposal(3)<N && proposal(4)<1)
        
        %log acceptance ratio (can be >1)
        lar = logLikelihoodSIR(proposal,data,N) - logLikelihoodSIR(theta,data,N);
        %prior ratio
        lar = lar + logPrior(proposal) - logPrior(theta);
        
        %generate uniform random number in (0,1)
        r = unifrnd(0,1);
        %accept if lar>log(r)
        if(lar>log(r))
            %proposal is new value of theta, increment accept counter
            theta = proposal;
            accept = accept+11;
            
        else
            %theta is unchanged, increment reject counter
            reject = reject+1;
        end
        
    %if not feasible    
    else
        %reject outside the range for theta
        reject = reject+1;
    end
    
    %store parameters for output every iteration
    stored(i,:) = theta;
    
    %recalculate and store trajectory for output
    para = struct("beta",theta(1),"gamma",theta(2),"iota",theta(3),"delta",theta(4),"N",N); 
    ICs = struct('S',para.N-theta(3),'I',theta(3),'R',0);
    Classes=ODE_SIR_model(para,ICs,max(data.t));
    incidence(i,:) = theta(4)*Classes.I/para.N;
    
    %update estimate of sigma using AM algorithm on multiples of n_0
    %i.e. every n_0 iterations, will only happen after n_0 iterations
    %(naturally)
    if (mod(i,n_0)==0)
        mu = mean(stored(1:i,:),1); %mean over 1st cpt
        sigma=cov(stored(1:i,:))+epsilon*eye(4);
    
    %this elseif only runs if the above if doesn't, so modulo doesn't need
    %to be checked
    elseif (i>n_0)
        muPrevious = mu;
        mu = (i*mu + stored(i,:)) / (i+1);
        sigma=(i-1)/i*sigma+(stored(i,:)'*stored(i,:)+i*muPrevious'*muPrevious-(i+1)*mu'*mu+epsilon*eye(4))/i;
        sigma = (sigma + sigma') / 2;
    end

end

%plots
figure(1)
clf

%Plot beta
subplot(4,1,1);
plot(1:iters,stored(:,1));
ylabel('$\beta$');
title(strcat('$\sigma$ = '," ",string(sigma(1,1)),', acceptance rate = ',string(accept/(accept+reject))));

%Plot gamma
subplot(4,1,2);
plot(1:iters,stored(:,2));
ylabel('$\gamma$');
title(strcat('$\sigma$ = '," ",string(sigma(2,2)),', acceptance rate = ',string(accept/(accept+reject))));

%Plot I(0)
subplot(4,1,3);
plot(1:iters,stored(:,3));
ylabel("I(0)");
title(strcat('$\sigma$ = '," ",string(sigma(3,3)),', acceptance rate = ',string(accept/(accept+reject))));

%plot delta
subplot(4,1,4);
plot(1:iters,stored(:,4));
ylabel("$\delta$");
title(strcat('$\sigma$ = '," ",string(sigma(4,4)),', acceptance rate = ',string(accept/(accept+reject))));


xlabel("iteration");


%Plot fit against data
figure(2)
clf
hold on
%First plot the shaded 95% area
x=0:max(data.t);
x2 = [x fliplr(x)];
y1=quantile(incidence,[0.975]); %upper CI
y3=quantile(incidence,[0.025]); %lower CI
Shaded = [y3 fliplr(y1)];
h2=fill(x2, Shaded, 'b','facealpha',0.3,'EdgeColor','none');

%Then plot median line
h1=plot(0:max(data.t),median(incidence,1),'b');

%Plot data as scatter
h3=scatter(data.t,data.x./data.n,'or');

legend([h1 h2 h3],'Median model incidence','95% CI model incidence','Data','Location','Southeast')
xlabel('Weeks')
ylabel('Incidence')
axis([0 max(data.t) 0 Inf])
hold off;

% return the whole output matrix plus accept and reject counters
output=struct('theta',stored,'incidence',incidence,'accept',accept,'reject',reject,'acceptanceRate',accept/(accept+reject));
    
    





