%Stochastic (tau leap) SIR model with demography code

function [Classes] = Tauleap_SEIRDmodel(para,ICs,maxtime,tau)


%Run the tauleap algorithm for an SIR model

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
V=ICs.V;
%Define the current time in the model as 0
t=0;

%Run the model until either the maxtime is exceeded or until there are no
%infected people remaining
while ((t<maxtime) && (E+I>0))

    %Define event rates 
    infection = para.beta*I*S/para.N;
    end_latency_period = para.sigma*E;
    recovery = (1-para.pd)*para.gamma*I;
    die = para.pd*para.gamma*I; %of measels
    %natural death
    Sdeath = para.mu*S;
    Edeath = para.mu*E;
    Ideath = para.mu*I;
    Rdeath = para.mu*R;
    Vdeath = para.mu*V;
    %births, 0.95 is pe
    immuneatbirth = 0.95*para.v*para.mu*para.N;
    notimmuneatbirth = ((1-para.v)+(1-0.95))*para.mu*para.N;
    
    
    %Compute how many events occur for each time step
    inf_events = poissrnd(tau*infection);
    lat_events = poissrnd(tau*end_latency_period);
    recovery_events = poissrnd(tau*recovery);
    die_events = poissrnd(tau*die);
    Sdeath_events = poissrnd(tau*Sdeath);
    Edeath_events = poissrnd(tau*Edeath);
    Ideath_events = poissrnd(tau*Ideath);
    Rdeath_events = poissrnd(tau*Rdeath);
    Vdeath_events = poissrnd(tau*Vdeath);
    immune_events = poissrnd(tau*immuneatbirth);
    nonimmune_events = poissrnd(tau*notimmuneatbirth);


    %Update events
    
    death_event_total = Sdeath_events + Edeath_events + Ideath_events + Rdeath_events + Vdeath_events;
    
    S = S + nonimmune_events - inf_events - Sdeath_events;
    E = E + inf_events - lat_events - Edeath_events;
    I = I + lat_events - recovery_events - die_events - Ideath_events;
    R = R + recovery_events - Rdeath_events;
    D = D + die_events;
    V = V + immune_events - Vdeath_events;


    %Check nothing less than zero
    %to be more rigorous, we need to distribute back events with correct
    %proportions because there are multiple possible events for each class
    %fortunately these cases are quite rare, so there isn't much extra
    %computation required
    if S<0 %"undo" exposure or natural death
        tmp=S;
        S=0;
        
        %to distribute events between the possible places they could have
        %gone to
        evtot = inf_events + Sdeath_events;
        
        %going to E
        ec = min(binornd(-tmp,inf_events/evtot),inf_events);
        %going back into S
        nc = -tmp-ec;
        
        %now put these back
        %undo exposure
        E = E + ec;
        %undo natural death
        S = S + nc;
    end
    if E<0 %"undo" exposure
        tmp=E;
        E=0;
        
        %to distribute events between the possible places they could have
        %gone to
        evtot = lat_events + Edeath_events;
        
        %going to I
        ic = min(binornd(-tmp,lat_events/evtot),lat_events);
        %going back into S
        nc = -tmp-ic;
        
        %now put these back
        %undo infection
        I = I + ic;
        %undo natural death
        S = S + nc;
    end
    if I<0 %"undo" recovery
        tmp=I;
        I=0;
        
        %to distribute events between the possible places they could have
        %gone to
        evtot = recovery_events+die_events+Ideath_events;
        
        %deaths to measles
        dc = min(binornd(-tmp,die_events/evtot),die_events);
        %recoveries
        rc = min(binornd(-tmp-dc,recovery_events/evtot),recovery_events);
        %natural deaths
        nc = -tmp-dc-rc;
        
        %now put these back
        %undo death
        D = D + dc;
        %undo recovery
        R = R + rc;
        %undo death/birth
        S = S + nc;
        
        %the commented out method is not that great
%         %I goes to R or D, so choose with rng based on their weighting
%         %if chosen class becomes negative, it is rectified later
%         if (rand(1)<para.pd && D>abs(tmp)) %D can't become less than 0
%             D=D+tmp; %take back out of D class
%         else
%             R=R+tmp; %take back out of R class
%         end
    end
    if R<0 %"undo" natural sdeath
        tmp=R;
        R=0;
        S=S+tmp; %take back out of S class
    end
    %D and V can't be less than 0, nobody can leave D or V classes


    %Update time
    t = t+tau;

    %Save information in the Classes structure by extending the matrices of the
    %model state and the associated time
    Classes.t = [Classes.t t];
    Classes.S = [Classes.S S];
    Classes.E = [Classes.E E];
    Classes.I = [Classes.I I];
    Classes.R = [Classes.R R];
    Classes.D = [Classes.D D];
    Classes.V = [Classes.V V];
    
end
    
