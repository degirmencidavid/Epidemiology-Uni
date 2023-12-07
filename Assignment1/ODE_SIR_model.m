%ODE SI model code

function [Classes] = ODE_SIR_model(para,ICs,maxtime)


%Run ODE using ODE45
opts = odeset("RelTol", 1e-5);
[t, pop] = ode45(@diff_SIRmodel,[0:1:maxtime],[ICs.S ICs.I ICs.R], opts, para);

%Convert output to structure
Classes = struct("S", pop(:,1),"I",pop(:,2),"R",pop(:,3),"t",t);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Diff equations

function dPop = diff_SIRmodel(t,pop,para)

%Assign the population matrix into the classes
S=pop(1);
I=pop(2);
R=pop(3);

%Write down the ODE system
dS = -para.beta*S*I/para.N;
dI = para.beta*S*I/para.N - para.gamma*I;
dR = para.gamma*I;

%Reshape the derivatives into a column vector
dPop = [dS; dI; dR];

end

end
