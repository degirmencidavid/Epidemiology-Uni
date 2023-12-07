%ODE SI model code

function [Classes] = ODE_SIR3_model(para,ICs,mintime,maxtime)


%Run ODE using ODE45
opts = odeset("RelTol", 1e-5);
[t, pop] = ode45(@diff_SIRmodel,[mintime:1:maxtime],[ICs.S ICs.Ir ICs.In ICs.R], opts, para);

%Convert output to structure
Classes = struct("S", pop(:,1),"Ir",pop(:,2),"In",pop(:,3),"R",pop(:,4),"t",t);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Diff equations

function dPop = diff_SIRmodel(t,pop,para)

%Assign the population matrix into the classes
S=pop(1);
Ir=pop(2);
In=pop(3);
R=pop(4);

%Write down the ODE system
dS = -para.beta*S*(Ir+In)/para.N;
dIr = para.p*para.beta*S*(Ir+In)/para.N - para.gamma*Ir;
dIn = (1-para.p)*para.beta*S*(Ir+In)/para.N - para.gamma*In;
dR = para.gamma*(Ir+In);

%Reshape the derivatives into a column vector
dPop = [dS; dIr; dIn; dR];

end

end
