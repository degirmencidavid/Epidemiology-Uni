%Stochastic (Gillespie) SIR model code

function [Classes] = Gillespie_SEIRDmodel(para,ICs,maxtime)


%Run Gillespie's algorithm for an SIR model

%Store the starting point of the simultion from the ICs and copy to a new
%structure called Classes. Define the starting time to be 0
Classes = ICs;
Classes.t = 0;

%Define the current state of the model from the ICs
S=ICs.S;
E=ICs.E;
I=ICs.I;
R=ICs.R;
D=ICs.D;
%Define the current time in the model as 0
t=0;

%Run the model until either the maxtime is exceeded or until there are no
%infected people remaining
while ((t<maxtime) && (E+I>0))

    %Define event rates 
    new_infection = para.beta*I*S/para.N; %i.e. exposed
    end_latency_period = para.sigma*E;
    recovery = (1-para.pd)*para.gamma*I;
    die = para.pd*para.gamma*I;

    %Compute the total rate
    R_total = new_infection + end_latency_period + recovery + die;

    %Compute time to next event using a uniform random number
    r_1=rand(1);
    T_nexttime = -log(r_1)/R_total;

    %Select which event has occured
    %First infection
    r_2=rand(1);

    if (r_2<new_infection/R_total)
        %Update infections
        E=E+1;
        S=S-1;
    elseif (r_2<(end_latency_period+new_infection)/R_total)
        %Update end latency
        I=I+1;
        E=E-1;
    elseif (r_2<(recovery+end_latency_period+new_infection)/R_total)
        %Update recoveries
        I=I-1;
        R=R+1;
    else
    %or it could be, more rigorously, but with the same effect:
    %elseif (r_2<(die+recovery+end_latency_period+new_infection)/R_total)
        %Update deaths
        I=I-1;
        D=D+1;
        %dead people removed from population
        %para.N = para.N - 1;
    end


    %Update time
    t = t+T_nexttime;

    %Save information in the Classes structure by extending the matrices of the
    %model state and the associated time
    Classes.t = [Classes.t t];
    Classes.S = [Classes.S S];
    Classes.E = [Classes.E E];
    Classes.I = [Classes.I I];
    Classes.R = [Classes.R R];
    Classes.D = [Classes.D D];
    
end