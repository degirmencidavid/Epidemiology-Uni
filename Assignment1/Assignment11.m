%Q1

%init

%Define model parameters as a structure
para = struct("beta",0.62,"gamma",1/2.6,"N",62000000);
%Define initial conditions as a structure
ICs = struct("S",para.N-1,"I",1,"R",0);
%Define time to run model for
maxtime = 250;
%Run model by calling function ODE_SIR_model.m
[Classes] = ODE_SIR_model(para,ICs,maxtime);

%(a) plot
figure(1);
plot(Classes.t, Classes.S, "color", [34 139 34]./255);
hold on
plot(Classes.t, Classes.I, "r", Classes.t, Classes.R, "b");
hold off
xlabel("Time / days");
ylabel("Number");
legend("S","I","R");
set(gca,"fontsize",16);
set(0,"defaultaxesfontsize",16);
set(0,"defaultlinelinewidth",2);

%(b)
%(i)The final size of epidemic Rinf
% Compute the final size of the outbreak using R(end) to the nearest person
Final_size = round(Classes.R(end));
% %Compute the final size of the outbreak using implict equation
%Final_size_eq = round(fsolve(@ (Rinf) [(para.N-1).*exp(-(R0/para.N)*Rinf) - para.N + Rinf], para.N/2));

%(ii)R0
R0 = para.beta/para.gamma;

%(iii)The epected outbreak duration
LastInft = find(Classes.I>1,1,"last");
Duration = ceil(Classes.t(LastInft));

%(c) Re(t)
FOI = para.beta.*Classes.S.*Classes.I./para.N;
Rnumber = FOI./(Classes.I*para.gamma);
figure(2)
plot(Classes.t,Rnumber);
xlabel("Time / days");
ylabel("R_e(t), reproduction number")

%(d)r(t)
rt = para.beta.*Classes.S.*Classes.I./para.N - para.gamma.*Classes.I;
figure(3)
plot(Classes.t,rt);
xlabel("Time / days");
ylabel("Growth rate, r(t)");

%alternative method
rt2 = para.gamma.*Classes.I.*(Rnumber - 1);
figure(4)
plot(Classes.t,rt2);
xlabel("Time / days");
ylabel("Growth rate, r(t)");
