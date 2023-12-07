%%%PREAMBLE%%%

clc
clear
clf

%addpath("..");

%formatting

%set default stuff, pretty sure half of this is redundant but w/e
set(0,'defaulttextinterpreter','latex')
set(0,'DefaultTextFontname', 'latex')
set(0,'DefaultAxesFontName', 'latex')
set(0,'DefaultLegendFontName', 'latex')
set(gca,"fontsize",16);
set(0,"defaultaxesfontsize",16);
set(0,"defaultlinelinewidth",2);


%%%END OF PREAMBLE%%%


%%READ CSV%%

%read "RSVdata.csv"
csvData = readtable("RSVdata.csv");
%convert to matrix
csvData = csvData{:,:};
%store columns in arrays
week = transpose(csvData(:,1));
hospitalisations = transpose(csvData(:,2));

%%END OF READ CSV%%



%%%INITIAL SET UP%%%

%population
N = 5500;

%starting guess for theta_0, these are per 10000 people
beta = 0.5;
gamma = 0.1;
iota = 3;
delta = 0.08;

%put these into theta_0
theta_0 = [beta gamma iota delta];
theta = theta_0;

%set up data as struct
Sample_Size = N;
n = Sample_Size*ones(1,length(hospitalisations));
data = struct('n',n,'t',week,'x',hospitalisations);
%%%END OF INITIAL SET UP%%%



%%%RUNNING MCMC ALGORITHM$$$
%run mcmc and put into myOutput
myOutput=mcmc(10000,theta,diag(theta)*10e-9,250,data);
%%%END OF RUNNING MCMC ALGORITHM



%%%HISTOGRAMS OF THETA%%%

figure(5)

%beta
subplot(2,2,1);
histogram(myOutput.theta(:,1));
xlabel("$\beta$");
ylabel("frequency");

%gamma
subplot(2,2,2);
histogram(myOutput.theta(:,2));
xlabel("$\gamma$");
ylabel("frequency");

%iota
subplot(2,2,3);
histogram(myOutput.theta(:,3));
xlabel("$\iota$");
ylabel("frequency");

%delta
%gamma
subplot(2,2,4);
histogram(myOutput.theta(:,4));
xlabel("$\delta$");
ylabel("frequency");

%%%END OF HISTOGRAMS OF THETA%%%



%%%ESTIMATE OF FINAL SIZE%%%

%theta values mean over 20 runs
thetam = [0.457 0.096 0.063 0.81];

%Define model parameters as a structure
para = struct("beta",thetam(1),"gamma",thetam(2),"N",5500);
%Define initial conditions as a structure
ICs = struct("S",para.N-theta(3),"I",theta(3),"R",0);
%Define time to run model for
maxtime = 47*7;
%Run model by calling function ODE_SIR_model.m
[Classes] = ODE_SIR_model(para,ICs,maxtime);

Final_size = round(Classes.R(end)*thetam(4));
%%%END OF ESTIMATE OF FINAL SIZE%%%


%%%NOTES ON ANONYMOUS FUNCTION METHOD TO CHOOSE STARTING POINT%%%
%anonymous function methods
%func = @(theta) log likelihood

%x = fminunc(logLikelihoodSIR(theta,data,N),theta)
%%%END OF NOTES%%%




%%%FOR FUN%%%

%fitting a polynomial to data (order 9)
figure(22)
p = polyfit(data.t,data.x,9);
x7 = linspace(1,10);
pl = polyval(p,data.t);
plot(data.t,pl)
hold on
scatter(data.t,data.x,'or');
hold off


%plot data, as proportion of total sample size
figure(23)
clf
scatter(data.t,data.x./data.n,'or');
hold on
plot(data.t,data.x./data.n,"--k");
hold on

%plot an estimation of poisson pdf, for fun
lambda = 32;
xl = data.t;
poiss = poisspdf(xl,lambda);
plot(xl,poiss);
hold off

%experimenting with betapdf
figure(24)
plot(0:0.01:1,betapdf(0:0.01:1,2,22))

%%%FUNTIME OVER%%%
