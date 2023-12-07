% Q1d - function to perform MCMC for all of the parameters
function [output] = mcmc2(iters,startingValue,sigma,data)

% define acceptance rate monitors
accept=0;
reject=0;

% current value for the MCMC is theta
theta=startingValue;
N=10000;

% prepare storage of output
stored=zeros(iters,length(theta));
prevalence=zeros(iters,1+max(data.t));

for i=1:iters
    % Generate a new proposal for beta
    proposal=mvnrnd(theta,sigma);
    % check proposal is in range
    if(min(proposal)>0 && proposal(3)<N && proposal(4)<1)
        % calculate log acceptance ratio (may be bigger than 1)
        lar = logLikelihoodSIR(proposal,data,N) - logLikelihoodSIR(theta,data,N);
        % don't forget the prior ratio!
        lar = lar + logPrior(proposal) - logPrior(theta);
        % generate a random number between 0 and 1;
        u = unifrnd(0,1);
        % accept if lar>log(u) iff ap>u 
        if lar>log(u)
            % the proposal becomes the new value of theta
            theta=proposal;
            accept=accept+1;
        else
            reject=reject+1;
        end
    else
        % automatically reject outside the range for beta (as it has prior
        % density zero).
        reject=reject+1;
    end
    % store parameters for output every iteration (at the moment)
    stored(i,:)=theta;
    % recalculate and store trajectory for output
    % (ok so I know I could do better/faster coding here, but I wanted to
    % keep the MH part of the MCMC as simple as possible).
    para = struct("beta",theta(1),"gamma",theta(2),"iota",theta(3),"delta",theta(4),"N",N); 
    ICs = struct('S',para.N-theta(3),'I',theta(3),'R',0);
    Classes=ODE_SIR_model(para,ICs,max(data.t));
    prevalence(i,:) = Classes.I/para.N; 
end


%plots
figure(3)
clf

%Plot beta
subplot(4,1,1);
plot(1:iters,stored(:,1));
ylabel('\beta');
title(strcat('\sigma = ',string(sigma(1,1)),', acceptance rate = ',string(accept/(accept+reject))));

%Plot gamma
subplot(4,1,2);
plot(1:iters,stored(:,2));
ylabel('\gamma');
title(strcat('\sigma = ',string(sigma(2,2)),', acceptance rate = ',string(accept/(accept+reject))));

%Plot I(0)
subplot(4,1,3);
plot(1:iters,stored(:,3));
ylabel("I(0)");
title(strcat('\sigma = ',string(sigma(3,3)),', acceptance rate = ',string(accept/(accept+reject))));

%plot delta
subplot(4,1,4);
plot(1:iters,stored(:,4));
ylabel("\delta");
title(strcat('\sigma = ',string(sigma(4,4)),', acceptance rate = ',string(accept/(accept+reject))));


xlabel("iteration");
%Plot fit against data
figure(4)
clf
hold on
%First plot the shaded 95% area
x=0:max(data.t);
x2 = [x fliplr(x)];
y1=0.9*quantile(prevalence,[0.975]); %upper CI
y3=0.9*quantile(prevalence,[0.025]); %lower CI
Shaded = [y3 fliplr(y1)];
h2=fill(x2, Shaded, 'b','facealpha',0.3,'EdgeColor','none');

%Then plot median line
h1=plot(0:max(data.t),0.9*median(prevalence,1),'b');

%Plot data as scatter
h3=scatter(data.t,data.x./data.n,'or');

legend([h1 h2 h3],'Median model prevalence','95% CI model prevalence','Data','Location','Southeast')
xlabel('Days')
ylabel('Prevalence')
axis([0 max(data.t) 0 Inf])
hold off;

% return the whole output matrix plus accept and reject counters
output=struct('theta',stored,'prevalence',prevalence,'accept',accept,'reject',reject,'acceptanceRate',accept/(accept+reject));
    
    

