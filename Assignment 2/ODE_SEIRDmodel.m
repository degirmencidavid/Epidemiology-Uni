%ODE SIR model code

function [Classes] = ODE_SEIRDmodel(para,ICs,maxtime)


%Run ODE using ODE45 using 1 day time steps
opts = odeset('RelTol',1e-5);
[t, pop] = ode45(@diff_SEIRDmodel, [0: 1: maxtime], [ICs.S ICs.E ICs.I ICs.R ICs.D], opts, para);

%Convert output to structure
Classes = struct('S',pop(:,1),'E',pop(:,2),'I',pop(:,3),'R',pop(:,4),'D',pop(:,5),'t',t);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Diff equations

    function dPop = diff_SEIRDmodel(t,pop,para)

        %Assign the population matrix into the classes
        S=pop(1);
        E=pop(2);
        I=pop(3);
        R=pop(4);
        D=pop(5);

        %Write down the ODE system
        dS = -para.beta*S*I/para.N;
        dE = para.beta*S*I/para.N - para.sigma*E;
        dI = para.sigma*E - para.gamma*I;
        dR = (1-para.pd) * para.gamma*I;
        dD = para.pd*para.gamma*I;
        
        %para.N = para.N - dD;??
        
        %Reshape the derivatives into a column vector
        dPop = [dS; dE; dI; dR; dD];

    end

end
