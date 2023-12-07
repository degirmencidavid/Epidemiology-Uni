%prior information we know/assume:
%R_0=beta/gamma
%delta is the proportion of infected individuals who are hospitalised,
%=>delta is in [0,1]

%From online resources:
%RSV R_0 is 2.4 to 3.6
%RSV infectious for 3-8 days, beta is in [1/8,1/3]
%RSV recovery time is 1-2 weeks, 7-14 days, gamma is in [1/14,1/7]
%RSV ~5.2/1000 infected children are hospitalised (delta)

%midpoints: R_0: 3.0, infectious period: 5.5 days, recovery time: 10.5 days
%beta = 1/5.5, gamma = 1/10.5, choose a different beta and gamma to better
%match up with R_0?

%delta: take mean to be 0.052, if no knowledge of this, uniform
%distribution between 0 and 1, otherwise distribution with mean 0.052

function [lp] = logPrior(theta)
    
    %gamma
    a_gamma = 10.5;
    b_gamma = 1;
    %lp = log(gampdf(theta(1), a_beta, b_beta)); %alternative
    %inverse gamma distribution method so we can use R0 information
    %to get gamma prior
    lp = -1/(b_gamma*theta(2))-(a_gamma+1)*log(theta(2));
    
    %alternative
    %%standard (non inverse) method for gamma prior
    %%gamma
    %a_beta = 5.5;
    %b_beta = 1;
    %%lp = lp + log(gampdf(theta(2), a_beta, b_beta));
    %%inverse \/
    %lp = lp -1/(b_gamma*theta(2))-(a_gamma+1)*log(theta(2));
    
    
    %R0 = beta/gamma
    a_r = 30;
    b_r = 1/10;
    lp = lp + log(gampdf(theta(1),a_r,b_r*theta(2)));
    
    
    %iota is unknown, choose a distribution at some point between 0 or 1 and
    %multiply by maximum population, we can't really know, so choose uniform
    %N*Beta(1,1) which is just a uniform distribution, 
    
    %alternative, but it's just 0 when you log it    
    %lp = lp + log(betapdf(theta(3),1,1));
    %log(unifpdf(theta(3),0,5500));
    lp = lp + 0;
   

    %delta is between 0 and 1
    %choose a distribution which is within that range - beta
    
    %has most of mass where we want it, not too informative either
    lp = lp + log(betapdf(theta(4),2,22));
    