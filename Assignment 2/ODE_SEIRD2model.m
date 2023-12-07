%ODE SIR model code

function [Classes] = ODE_SEIRD2model(para,ICs,maxtime)


%Run ODE using ODE45 using 1 day time steps
opts = odeset('RelTol',1e-5);
[t, pop] = ode45(@diff_SIRmodel, [0: 1: maxtime], [ICs.S ICs.E ICs.I ICs.R ICs.D ICs.V], opts, para);

%Convert output to structure
Classes = struct('S',pop(:,1),'E',pop(:,2),'I',pop(:,3),'R',pop(:,4),'D',pop(:,5),'V',pop(:,6),'t',t);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Diff equations

    function dPop = diff_SIRmodel(t,pop,para)

        %Assign the population matrix into the classes
        S=pop(1);
        E=pop(2);
        I=pop(3);
        R=pop(4);
        D=pop(5);
        V=pop(6);

        %Write down the ODE system
        dS = ((1-0.95)*para.v+(1-para.v))*para.mu*para.N - para.beta*S*I/para.N - para.mu*S;
        dE = para.beta*S*I/para.N - para.sigma*E - para.mu*E;
        dI = para.sigma*E - para.gamma*I - para.mu*I;
        dR = (1-para.pd) * para.gamma*I - para.mu*R;
        dD = para.pd*para.gamma*I - para.mu*D;
        dV = 0.95*para.v*para.mu*para.N - para.mu*V;
        
        %para.N = (S+E+I+R);?????

        %Reshape the derivatives into a column vector
        dPop = [dS; dE; dI; dR; dD; dV];

    end

end
