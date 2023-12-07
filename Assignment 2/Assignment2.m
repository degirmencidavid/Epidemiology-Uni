%%%PREAMBLE%%%

%clear workspace and cmd window
clear
clc


%set default stuff, pretty sure half of this is redundant but w/e
set(0,'defaulttextinterpreter','latex')
set(0,'DefaultTextFontname', 'latex')
set(0,'DefaultAxesFontName', 'latex')
set(0,'DefaultLegendFontName', 'latex')
set(gca,"fontsize",16);
set(0,"defaultaxesfontsize",16);
set(0,"defaultlinelinewidth",2);

%%%END OF PREAMBLE%%%




%%%QUESTION 1%%%

%immune proportion
immune = 0.9;

%%Question 1 (a)%%

%initialisation of model

%Define model parameters as a structure
para = struct("beta",5,"pd",0.01,"sigma",1/10,"gamma",1/3,"N",100000); 
%Define initial conditions as a structure
ICs = struct("S",para.N-1,"E",0,"I",1,"R",0,"D",0);
%Define time to run model for (2 years like in Q2)
maxtime = 730;
%Run model by calling function ODE_SIRmodel.m
[Classes] = ODE_SEIRDmodel(para,ICs,maxtime);


%%with no immunity%%

%R0 is unchanged - no demography, R0 = beta/gamma
R0 = para.beta/para.gamma;

%Re(0) = beta*S(0)/gamma*N, S(0) = N, so it follows Re(0) = R0
Re0 = R0;

%plot of Re for my sanity%

%Re(t) = beta*S(t)/gamma*N
Re = para.beta.*Classes.S./(para.gamma*para.N);

%the plot itself
figure(1)
clf
plot(Classes.t, Re);
xlabel("Time / days");
ylabel("R${_e}$(t)");

%end of Re plot%

%%end of with no immunity%%




%%with immunity%%


%reinitialise model with immunity being considered

%define these new structs and variables with I appended to the end
%to represent Immunity


%immune proportion
immune = 0.9;

%Define model parameters as a structure
paraI = struct("beta",5,"pd",0.01,"sigma",1/10,"gamma",1/3,"N",100000); 
%Define initial conditions as a structure
ICsI = struct("S",(1-immune)*paraI.N-1,"E",0,"I",1,"R",immune*paraI.N,"D",0);
%Define time to run model for (2 years)
maxtimeI = 730;
%Run model by calling function ODE_SIRmodel.m
[ClassesI] = ODE_SEIRDmodel(paraI,ICsI,maxtimeI);


%calculating R0 with immunity
%S(0) * contact rate * avg time spent infectious
%=(1-immune)*N * (beta/N) * 1/gamma
%=(1-immune)*beta/gamma
RI0 = (1-immune)*paraI.beta/paraI.gamma;


%Re(t) = beta*S(t)/gamma*(1-immune)*N
ReI = paraI.beta.*ClassesI.S/(paraI.gamma*paraI.N);

%but, as Re(0)=R0
%Re(1) of immune case
ReI1 = ReI(1);
%Re(0) of immune case
ReI0 = RI0;

%plot of Re(t) with immunity%
figure(2)
clf
plot(ClassesI.t, ReI);
%end of plot of Re(t) with immunity%

%this plot shows the behaviour of Re(t) with immunity, we can see that
%as t->0, Re(t)->1.5 which is what I calculated before

%%end of with immunity%%




%(b)

%initial conditions ((1-immune)*N,0,1,0,0)

%plot SEIRD with immunity
figure(3);
clf
plot(ClassesI.t, ClassesI.S, "color", [34 139 34]./255);
hold on
plot(ClassesI.t, ClassesI.E, "m",ClassesI.t, ClassesI.I,"r")
hold on
plot(ClassesI.t, ClassesI.R, "b", ClassesI.t, ClassesI.D,"c");
hold off
xlabel("Time / days");
ylabel("Number");
legend("S","E","I","R","D",'Interpreter','latex');



%plot I individually
figure(4)
clf
plot(ClassesI.t, ClassesI.I, "r");
xlabel("Time / days");
ylabel("Number of Infections (I)");
legend("I",'Interpreter','latex');



%total amount of infected individuals
Itot = trapz(ClassesI.t, ClassesI.I);
%deaths are cumulative
Dtot = paraI.pd*paraI.gamma*Itot;
Dtot_by_end = ClassesI.D(end);
%tiny difference of 0.0001 - likely due to float miscalculation

%doesn't matter which Dtot I use, rounding will give 58 either way
rounded_deaths = round(Dtot_by_end);

%it is not this:
%DtotCalc = trapz(ClassesI.t, ClassesI.D);
%because D is cumulative


%plot D individually
figure(5)
clf
plot(ClassesI.t, ClassesI.D);
xlabel("Time / days");
ylabel("Cumulative Number of Deaths (D)");
legend("D",'Interpreter','latex');



%%%END OF QUESTION 1%%%






%%%QUESTION 2%%%

%2(b)

%Run stochastic model using Gillespie once
[ClassesG] = Gillespie_SEIRDmodel(paraI,ICsI,maxtimeI);

%plot single gillespie run on figure with ode
figure(6)
clf
plot(ClassesG.t,ClassesG.I,"r",ClassesI.t,ClassesI.I,"--k")
axis([0 730 0 Inf])
xlabel("Time / days");
ylabel("Infections (I)");
legend("Gillespie","ODE");

%zoom in to show this better (only if no take off because it is the most
%likely case), this can be set appropriately if needed
% axis([0 10 0 3])


%(c)

%some properties of Gillespie model

%Peak size
Peak_height_I=max(ClassesG.I);
%Final size
Final_size_I = ClassesG.R(end)+ClassesG.D(end)-immune*paraI.N.;
%Duration
[,ix]=find(ClassesG.I+ClassesG.E>=1,1,"last");
Duration_I=ClassesG.t(ix);


%(d)

immune = 0.9;
%%plotting multiple gillespie runs on one figure with ode to compare%%

%Define model parameters as a structure
paraI = struct("beta",5,"pd",0.01,"sigma",1/10,"gamma",1/3,"N",100000); 
%Define initial conditions as a structure
ICsI = struct("S",(1-immune)*paraI.N-1,"E",0,"I",1,"R",immune*paraI.N,"D",0);
%Define time to run model for (2 years)
maxtimeI = 730;

%set up plot 
figure(7)
clf
%plot(ClassesG.t,ClassesG.I,"r",ClassesI.t,ClassesI.I,"--k")
hold on

%Use a for loop to run each simulation, store metrics then plot the
%dynamics
NRuns=10;
tic
for i=1:NRuns
    [Classes_G] = Gillespie_SEIRDmodel(paraI,ICsI,maxtimeI);
    
    %Store stochastic simulation metrics

    %Peak size
    Peak_height(i)=max(Classes_G.I);
    %Final size is recovered + deaths
    Final_size(i) = Classes_G.R(end)+Classes_G.D(end)-immune*paraI.N;
    %Duration
    [,ix]=find(Classes_G.I+Classes_G.E>=1,1,"last");
    Duration(i)=Classes_G.t(ix);
    
    h1=plot(Classes_G.t,Classes_G.I,"r");
end
toc
%Add the ODE solution back on top of the stochastic solutions
h2=plot(ClassesI.t,ClassesI.I,"--k");
hold off
axis([0 730 0 Inf])
xlabel("Time / days");
ylabel("Infections (I)");
legend([h1 h2], "Gillespie","ODE");

%%end of plot%%


%d needs some work

%%Some ode properties%%

%Peak size
Peak_height_ODE=max(ClassesI.I);
%Final size is recovered + deaths
Final_size_ODE = ClassesI.R(end)+ClassesI.D(end)-immune*paraI.N;
%Duration (remember we have to define a proxy for when the epidemic ends)
[,ix]=find(ClassesI.I>1,1,"last");
Duration_ODE=ClassesI.t(max([ix 1]));

%%end of ode properties%%




%%%PLOT OF HISTOGRAMS%%%

%Plot histograms of: Size of peak, duration and final size

%Size of peak
figure(8)
clf
histogram(Peak_height,"numbins",30,"normalization","probability");
xlabel("Peak height");
ylabel("Probability");
hold on
ymax = get(gca,"ylim");
plot([Peak_height_ODE Peak_height_ODE], [0 ymax(2)]);
hold off
legend("1000 Gillespie realisations","ODE");

%Size of outbreak
figure(9)
clf
histogram(Final_size,"numbins",50,"normalization","probability");
xlabel("Final size");
ylabel("Probability");
hold on
ymax = get(gca,"ylim");
plot([Final_size_ODE Final_size_ODE], [0 ymax(2)]);
hold off
legend("1000 Gillespie realisations","ODE");
%can make it look nicer with edges rather than numbins

%Duration
figure(10)
clf
histogram(Duration,"numbins",30,"normalization","probability");
xlabel("Duration / days");
ylabel("Probability");
hold on
ymax = get(gca,"ylim");
plot([Duration_ODE Duration_ODE], [0 ymax(2)]);
hold off
legend("1000 Gillespie realisations","ODE");


%%%END OF HISTOGRAM PLOTS%%%



%(e)


%%analytically%%

%to become infected, one must first become exposed

%total of possible cases
Tot = paraI.beta+paraI.gamma;

%probability infected
probI = (paraI.beta)/Tot;

%probability recover
probR = paraI.gamma/Tot;

%0 extra infections, P(I)^0 * P(R)^1 * number of possible combinations
%of paths to get to 0 infections, which is 1, so just P(R)
prob0 = probR;

%1 extra infection, again only 1 possible path, infect->recover->recover
%P(I)^1 * P(R)^2 * 1
prob1 = probI*probR.^2;

%2 extra infections,
%possible combinations are:
%infect->infect->recover->recover->recover
%infect->recover->infect->recover->recover
%P(I)^2 * P(R)^3 * 2
prob2 = probI.^2*probR.^3 * 2;

%3 extra infections
%possible combinations are:
%infect->recover->infect->recover->infect->recover->recover
%infect->infect->recover->infect->recover->recover->recover
%infect->recover->infect->infect->recover->recover->recover
%infect->infect->infect->recover->recover->recover->recover
%infect->infect->recover->recover->infect->recover->recover
%5 of them, so P(I)^3 * P(R)^4 * 5
prob3 = probI.^3*probR.^4 * 5;

%probability of between 0-3 extra infections = sum of these
probTot = prob0 + prob1 + prob2 + prob3;

%%end of analytically%%


%%via model%%

%final size of outbreak is between 1 and 4
%total number of "outbreaks" that only have 0-3 additional infections
sum_outbreak = 0;
for i=1:NRuns
    %check final size is between 1 and 4 and add 1 to the sum
    %the check for Final_size(i)>=1 isn't actually necessary
    %but it remains to be more understandable
    if (Final_size(i)<=4 && Final_size(i)>=1)
        sum_outbreak = sum_outbreak + 1;
    end
end

%proportion that there are between 1 and 4 final size
proportion14 = sum_outbreak/NRuns;

%%end of via model%%



%%%END OF QUESTION 2%%%





%%%QUESTION 3%%%


%(b)

%tau leap

%defining some variables

%proportion of infants vaccinated, 0% here
vaxrate = 0;
%1 million population
population = 1000000;
%1 day timestep
timestepTL = 1;

%Define model parameters as a structure
paraTL = struct("beta",5,"pd",0.01,"sigma",1/10,"gamma",1/3,"mu",1/(60*365),"v",vaxrate,"N",population); 
%Define initial conditions as a structure
ICsTL = struct("S",66700,"E",430,"I",130,"R",751270,"D",0,"V",181470);
%Define time to run model for, 20 years
maxtimeTL = 7300;
%Run model by calling function ODE_SIRmodel.m
[ClassesTL] = Tauleap_SEIRDmodel(paraTL,ICsTL,maxtimeTL,timestepTL);

%%%TEST%%%
% testing for original conditions
% paraTL = struct("beta",5,"pd",0.01,"sigma",1/10,"gamma",1/3,"mu",0,"v",0,"N",100000);
% ICsTL = struct("S",0.1*paraTL.N-1,"E",0,"I",1,"R",0,"D",0,"V",0);
% [ClassesTL] = Tauleap_SEIRDmodel(paraTL,ICsTL,maxtimeTL,timestepTL);
%%%END OF TEST%%%


%ODE 
%Run model by calling function ODE_SEIRD2model.m
[ClassesODE] = ODE_SEIRD2model(paraTL,ICsTL,maxtimeTL);

%plot single tau leap
figure(11)
clf
plot(ClassesTL.t./365,ClassesTL.I,"r")
hold on
plot(ClassesODE.t./365, ClassesODE.I,"--k");
hold off
xlabel("Time / years");
ylabel("Infections (I)");
legend("Tau leap","ODE");


%ODE 
%Run model by calling function ODE_SEIRD2model.m
[ClassesODE] = ODE_SEIRD2model(paraTL,ICsTL,maxtimeTL);

%plotting many tau leap runs with ode
figure(12)
clf
hold on

%Use a for loop to run each simulation
%dynamics
NRuns=10;
%intialise empty arrays - better memory management
deathruns = zeros(1,NRuns);
durationstl = zeros(1,NRuns);
tic
for i=1:NRuns
    
    %calculate each run
    [Classes_TL] = Tauleap_SEIRDmodel(paraTL,ICsTL,maxtimeTL,timestepTL);
        
    %measles deaths
    deathruns(i) = Classes_TL.D(end);
    
    %duration
    [,ix]=find(Classes_TL.I>1,1,"last");
    Duration_TL_1=Classes_TL.t(max([ix 1]));
    durationstl(i) = Duration_TL_1;
       
    %plot each run
    h1=plot(Classes_TL.t./365,Classes_TL.I,'color','r');
   
end
toc

%%%TEST%%%
% %Define model parameters as a structure
% paraI = struct("beta",5,"pd",0.01,"sigma",1/10,"gamma",1/3,"N",100000); 
% %Define initial conditions as a structure
% ICsI = struct("S",(1-immune)*paraI.N-1,"E",0,"I",1,"R",0,"D",0);
% %Define time to run model for (2 years)
% maxtimeI = 730;
% %Run model by calling function ODE_SIRmodel.m
% [ClassesI] = ODE_SEIRDmodel(paraI,ICsI,maxtimeI);
% %Add the ODE solution back on top of the stochastic solutions
% plot(ClassesI.t/365,ClassesI.I,"--k");
%%%END OF TEST%%%



%plot deathruns as histogram
figure(22)
histogram(deathruns,"numbins",30,"normalization","probability");
xlabel("Number of deaths due to measles");
ylabel("Probability");

%plot durations as histogram
figure(23)
histogram(durationstl./365,"numbins",40,"normalization","probability");
xlabel("Duration / years");
ylabel("Probability");




%(d)

%again, setting up variables
%redundant for some, but it is for readability and so I can control
%stuff more easily to experiment with them
%proportion of infants vaccinated, 60% here
vaxrate = 0.6;
%1 million population
population = 1000000;
%1 day timestep
timestepTL = 1;

%Define model parameters as a structure
paraTL = struct("beta",5,"pd",0.01,"sigma",1/10,"gamma",1/3,"mu",1/(60*365),"v",vaxrate,"N",population); 
%Define initial conditions as a structure
ICsTL = struct("S",66700,"E",430,"I",130,"R",751270,"D",0,"V",181470);
%Define time to run model for, 20 years
maxtimeTL = 7300;

%Run TL model by calling function Tauleap_SEIRDmodel.m
[ClassesTL] = Tauleap_SEIRDmodel(paraTL,ICsTL,maxtimeTL,timestepTL);




%plot multiple tau leap simulations
figure(13)
clf
hold on

%Use a for loop to run each simulation, store metrics then plot the
%dynamics

%initial set up for probabilities
%sum for less than 5 years
lt5y = 0;
%sum for greater than 5 years
gt5y = 0;

%loop to run tau leap NRuns times - very long
NRuns=1000;
tic
for i=1:NRuns
    
    %calculate each run
    [Classes_TL] = Tauleap_SEIRDmodel(paraTL,ICsTL,maxtimeTL,timestepTL);
    
    %plot each run
    h1=plot(Classes_TL.t./365,Classes_TL.I,'color','r');
    
    %measles deaths
    deathruns(i) = Classes_TL.D(end);
    
    %duration
    [,ix]=find(Classes_TL.I>1,1,"last");
    Duration_TL_1=Classes_TL.t(max([ix 1]));
    durationstl(i) = Duration_TL_1;
    
    %set up for probability of elimination
    if (Classes_TL.t(end)/365<5)
        lt5y = lt5y+1;
    else
        gt5y = gt5y + 1;
    end
   
end
toc
hold on
%Run model by calling function ODE_SEIRD2model.m
[ClassesODE] = ODE_SEIRD2model(paraTL,ICsTL,maxtimeTL);

%plot ode
h2=plot(ClassesODE.t./365,ClassesODE.I,"--k");
xlabel('Time / years')
ylabel('Infections (I)')
legend([h1,h2],'Tau leap','ODE')
hold off

%plot deathruns as histogram
figure(24)
histogram(deathruns,"numbins",30,"normalization","probability");
xlabel("Number of deaths due to measles");
ylabel("Probability");

%plot durations as histogram
figure(25)
histogram(durationstl./365,"numbins",40,"normalization","probability");
xlabel("Duration / years");
ylabel("Probability");



%%probabilities

%probability  of elimination in less than 5 years
problt5y = lt5y/NRuns;

%probability  of elimination in greater than 5 years
probgt5y = gt5y/NRuns;



%%%CHECK%%%
%check if something's wrong with total of all classes
%plot ode sum of all classes, should still be 1 million, which it is
figure(14)
clf
plot(ClassesODE.t,ClassesODE.S+ClassesODE.E+ClassesODE.I+ClassesODE.R+ClassesODE.V+ClassesODE.D);
%%%END OF CHECK%%%


%to work out what vaccination coverage would be needed can be done
%a few ways, trial and error is the easiest

%sweet spot is 65 percent (to the nearest 1 percent)
paraTL.v = 0.84;
    
%initial set up for probabilities
%sum for less than 5 years
lt5y2 = 0;
%sum for greater than 5 years
gt5y2 = 0;
%loop to run tau leap NRuns times - very long
NRuns=1000;
for i=1:NRuns
    
    %calculate each run
    [Classes_TL] = Tauleap_SEIRDmodel(paraTL,ICsTL,maxtimeTL,timestepTL);
    
    %plot isn't needed
    
    %set up for probability of elimination
    if (Classes_TL.t(end)/365<5)
        lt5y2 = lt5y2+1;
    else
        gt5y2 = gt5y2 + 1;
    end
   
end


%%probabilities

%probability  of elimination in less than 5 years
problt5y2 = lt5y2/NRuns;

%probability  of elimination in greater than 5 years
probgt5y2 = gt5y2/NRuns;




