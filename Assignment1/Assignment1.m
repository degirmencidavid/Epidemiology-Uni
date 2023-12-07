%preamble, clear workspace and cmd window
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


%NB - para,ICs, etc are recreated each time for readability


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

%(a) plot SIR 
figure(1);
plot(Classes.t, Classes.S, "color", [34 139 34]./255);
hold on
plot(Classes.t, Classes.I, "r", Classes.t, Classes.R, "b");
hold off
xlabel("Time / days");
ylabel("Number");
legend("S","I","R",'Interpreter','latex');

%(b)
%(i)The final size of epidemic Rinf
% Compute the final size of the outbreak using R(end) to the nearest person
Final_size = round(Classes.R(end));
% %Compute the final size of the outbreak using implict equation
%Final_size_eq = round(fsolve(@ (Rinf) [(para.N-1).*exp(-(R0/para.N)*Rinf) - para.N + Rinf], para.N/2));

%(ii)R0
R0 = para.beta/para.gamma;

%(iii)The epected outbreak duration
LastInft = find(Classes.I>=1,1,"last");
Duration1 = ceil(Classes.t(LastInft));

%(c) Re(t)
FOI = para.beta.*Classes.S.*Classes.I./para.N;
Rnumber = FOI./(Classes.I*para.gamma);
figure(2)
plot(Classes.t,Rnumber);
xlabel("Time / days");
ylabel("R$_e$(t), reproduction number")

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


%Q2

%init

%Define model parameters as a structure
para2 = struct("beta",0.62,"p",0.04,"gamma",1/2.6,"N",62000000);
%Define initial conditions as a structure
ICs2 = struct("S",0.87*para2.N-1,"Ir",0,"In",1,"R",0);
%Define time to run model for
maxtime2 = 250;
%Run model by calling function ODE_SIR_model.m
[Classes2] = ODE_SIR2_model(para2,ICs2,maxtime2);

%plot sir
figure(5);
plot(Classes2.t, Classes2.S, "color", [34 139 34]./255);
hold on
plot(Classes2.t, Classes2.In + Classes2.Ir, "r", Classes2.t, Classes2.R, "b");
hold off
xlabel("Time / days");
ylabel("Number");
legend("S","I","R",'Interpreter','latex');


% %plot sirinr
% figure(8);
% plot(Classes2.t, Classes2.S, "color", [34 139 34]./255);
% hold on
% plot(Classes2.t, Classes2.In,"g",Classes2.t, Classes2.Ir, "r", Classes2.t, Classes2.R, "b");
% hold off
% xlabel("Time / days");
% ylabel("Number");
% legend("S","In","Ir","R");
% set(gca,"fontsize",16);
% set(0,"defaultaxesfontsize",16);
% set(0,"defaultlinelinewidth",2)

%force of infection
FOI2 = para2.beta.*(Classes2.Ir+Classes2.In).*Classes2.S./para2.N;

%cumulative total reported cases
trc = cumtrapz(Classes2.t,para2.p.*FOI2);

%daily incidence
dit(1) = trc(1);
for i=2:size(trc)
    dit(i) = trc(i) - trc(i-1);
end

%plot daily incidence
figure(6);
plot(Classes2.t,round(dit));
xlabel("Time / days");
ylabel("Daily Reported Incidence");

%plot cumulative number of reported infections
figure(7);
plot(Classes2.t,round(trc))
xlabel("Time / days");
ylabel("Total Reported Cases");

round(trc(end))

%check if cumulative reported cases makes sense (they do :), very nice)
%disp(round(trc(end)*25)-round(Classes2.R(end)))




%Q3

%a
%import csv
tablg = readtable("SwineFluCaseData.csv");
dates = tablg{:,2};
cases = tablg{:,4};
cumcases = tablg{:,5};

%bar plot of weekly reported cases
figure(9);
bar(dates-24,cases);
xlabel("Week");
ylabel("Weekly Cases");

figure(10);
bar(dates-24,cumcases);
xlabel("Week");
ylabel("Weekly Cases Cumulative");

%b

%

%infection starts at beginning of week 17 = day 119, beta changes on day
%210
mintime3 = 119;
maxtime3 = 210;

%Define model parameters as a structure
para3 = struct("beta",0.62,"p",0.04,"gamma",1/2.6,"N",62000000);
%Define initial conditions as a structure
ICs3 = struct("S",0.87*para3.N-1,"Ir",0,"In",1,"R",0);
%set up model
[Classes3] = ODE_SIR3_model(para3,ICs3,mintime3,maxtime3);

%summer starts, betaholidays = 0.6beta
mintime4 = 210;
maxtime4 = 245;

%Define model parameters as a structure
para4 = struct("beta",0.62*0.6,"p",0.04,"gamma",1/2.6,"N",62000000);
%new ics are last conditions from Classes3
ICs4 = struct("S",Classes3.S(end),"Ir",Classes3.Ir(end),"In",Classes3.In(end),"R",Classes3.R(end));
[Classes4] = ODE_SIR3_model(para4,ICs4,mintime4,maxtime4);

%after school holidays, beta returns to original
mintime5 = 245;
maxtime5 = 500;

para5 = struct("beta",0.62,"p",0.04,"gamma",1/2.6,"N",62000000);
%new ics
ICs5 = struct("S",Classes4.S(end),"Ir",Classes4.Ir(end),"In",Classes4.In(end),"R",Classes4.R(end));
%set up model
[Classes5] = ODE_SIR3_model(para5,ICs5,mintime5,maxtime5);

%plot
figure(11);
plot(Classes3.t-119,Classes3.Ir+Classes3.In);
hold on
plot(Classes4.t-119,Classes4.Ir+Classes4.In);
hold on
plot(Classes5.t-119,Classes5.Ir+Classes5.In);
hold off
xlabel("Time / days");
ylabel("Number of infections, I(t)=Ir(t)+In(t)");
l11 = legend("Pre-holiday, $\beta=0.62$","During holiday, $\beta=0.372$","Post-holiday, $\beta=0.62$",'Interpreter','latex');
set(l11, 'interpreter', 'latex')
%c
Re3 = para3.beta.*Classes3.S./(para3.gamma*para3.N);
Re4 = para4.beta.*Classes4.S./(para4.gamma*para4.N);
Re5 = para5.beta.*Classes5.S./(para5.gamma*para5.N);

%plot R_e(t)
figure(12);
plot(Classes3.t-119,Re3);
hold on
plot(Classes4.t-119,Re4);
hold on
plot(Classes5.t-119,Re5);
hold on
plot([Classes3.t(end)-119,Classes4.t(1)-119],[Re3(end),Re4(1)],"--");
hold on
plot([Classes4.t(end)-119,Classes5.t(1)-119],[Re4(end),Re5(1)],"--");
hold off
xlabel("Time / days");
ylabel("R$_e$(t), reproduction number","interpreter","latex");
l12 = legend("Pre-holiday, $\beta=0.62$","During holiday, $\beta=0.372$","Post-holiday, $\beta=0.62$",'Interpreter','latex');
set(l12, 'interpreter', 'latex')

%d

%final size by plot
fin_size = round(Classes5.R(end));

%Duration of outbreak, taken as when In+Ir<1
LastInft = find(Classes5.Ir+Classes5.In>=1,1,"last");
Duration = ceil(Classes5.t(LastInft))-119;



%e 1 method in example sheet 2

%calculate weekly incidence
FOI3 = para3.beta.*(Classes3.Ir).*Classes3.S./para3.N;
for it=1:floor(size(FOI3)/7)
    Weekly_inc3(it) = trapz(Classes3.t([7*(it-1):7*it]+1), FOI3([7*(it-1):7*it]+1));
end

FOI4 = para4.beta.*(Classes4.Ir).*Classes4.S./para4.N;
for it=1:floor(size(FOI4)/7)
    Weekly_inc4(it) = trapz(Classes4.t([7*(it-1):7*it]+1), FOI4([7*(it-1):7*it]+1));
end

FOI5 = para5.beta.*(Classes5.Ir).*Classes5.S./para5.N;
for it=1:floor(size(FOI5)/7)
    Weekly_inc5(it) = trapz(Classes5.t([7*(it-1):7*it]+1), FOI5([7*(it-1):7*it]+1));
end

wi = [Weekly_inc3,Weekly_inc4,Weekly_inc5];


%plot
figure(13);
bar(dates-17,cases);
hold on
plot([17:70]-17,wi);
hold off
xlabel("Time / weeks");
ylabel("Weekly reported incidence");
legend("Actual data", "Model Weekly incidence",'Interpreter','latex');



%e 2 by days scaled down, this one is wrong but I left it in
FOI3 = para3.beta.*(Classes3.Ir+Classes3.In).*Classes3.S./para3.N;
trc3 = cumtrapz(Classes3.t,para3.p.*FOI3);
%daily incidence
dit3(1) = trc3(1);
for i=2:size(trc3)
    dit3(i) = trc3(i) - trc3(i-1);
end


FOI4 = para4.beta.*(Classes4.Ir+Classes4.In).*Classes4.S./para4.N;
trc4 = cumtrapz(Classes4.t,para4.p.*FOI4);
%daily incidence
dit4(1) = dit3(end);
for i=2:size(trc4)
    dit4(i) = trc4(i) - trc4(i-1);
end

FOI5 = para5.beta.*(Classes5.Ir+Classes5.In).*Classes5.S./para5.N;
trc5 = cumtrapz(Classes5.t,para5.p.*FOI5);
%daily incidence
dit5(1) = dit4(end);
for i=2:size(trc5)
    dit5(i) = trc5(i) - trc5(i-1);
end

%connect to 1 list
days = transpose([dit3,dit4,dit5]);
times = [Classes3.t;Classes4.t;Classes5.t];

it=1;
% while(it<=size(days,1))
%     sum = 0;
%     if(it+6<=size(days,1))
%         for i=it:it+7
%             sum = sum + days(i);
%         end
%     end
%     wit(it) = sum;
%     it = it+7;
% end

%sum 7 days, then go to next week do the same etc.
while(it<size(days,1))
    sum  = 0;
    if(it+7<=size(days,1))
        for i=1:7
            sum = sum+days(it+i);
        end
    end
    it = it+7;
    wit(it) = sum;
end

%because I'm lazy, wit has elements every 7 for weekly incidence, so
%I can just remove the junk 0s like this
wit = nonzeros(wit);

%plot
figure(14);
bar(dates-17,cases);
hold on
plot([17:17+size(wit)-1]-17,wit);
hold off
xlabel("Time / weeks");
ylabel("Weekly reported incidence vs data collected");
legend("Actual data", "Model Weekly incidence",'Interpreter','latex');


%e 3
FOI3 = para3.beta.*(Classes3.Ir+Classes3.In).*Classes3.S./para3.N;
trc3 = cumtrapz(Classes3.t,para3.p.*FOI3);

var3 = size(trc3);
i=8;
j=1;
while(i<var3(1)+1)
   s = [var3(1),i];
   wdit3(j) =  trc3(i) - trc3(i-7);
   j = j+1;
   i = i+7;
end


FOI4 = para4.beta.*(Classes4.Ir+Classes4.In).*Classes4.S./para4.N;
trc4 = cumtrapz(Classes4.t,para4.p.*FOI4);

var4 = size(trc4);
i=8;
j=1;
while(i<var4(1)+1)
   s = [var4(1),i];
   wdit4(j) =  trc4(i) - trc4(i-7);
   j = j+1;
   i = i+7;
end


FOI5 = para5.beta.*(Classes5.Ir+Classes5.In).*Classes5.S./para5.N;
trc5 = cumtrapz(Classes5.t,para5.p.*FOI5);

var5 = size(trc5);
i=8;
j=1;
while(i<var5(1)+1)
   s = [var5(1),i];
   wdit5(j) =  trc5(i) - trc5(i-7);
   j = j+1;
   i = i+7;
end

ditot = [wdit3,wdit4,wdit5];
sizzle = size(ditot);

figure(15);
bar(dates-17,cases);
hold on
plot([17:70]-17,ditot);
hold off
xlabel("Time / weeks");
ylabel("Weekly reported incidence vs data collected");
legend("Actual data", "Model Weekly incidence",'Interpreter','latex');


%Q4

%case in Q3

%latex graphs here please


%paras and ics with various conditions
%with initial immunity
para6 = struct("beta",0.62,"p",0.04,"gamma",1/2.6,"N",62000000);
ICs6 = struct("S",0.87*para6.N-1,"Ir",0,"In",1,"R",0);
mintime6 = 119;
maxtime6 = 210;
[Classes6] = ODE_SIR3_model(para6,ICs6,mintime6,maxtime6);

R06 = para6.beta/para6.gamma;
HIT6perc = 100*(1-(1/R06));

para7 = struct("beta",0.62*0.6,"p",0.04,"gamma",1/2.6,"N",62000000);
ICs7 = struct("S",Classes6.S(end),"Ir",Classes6.Ir(end),"In",Classes6.In(end),"R",Classes6.R(end));
mintime7 = 210;
maxtime7 = 245;
[Classes7] = ODE_SIR3_model(para7,ICs7,mintime7,maxtime7);

R07 = para7.beta/para7.gamma;
HIT7perc = 100*(1-(1/R07));

para8 = struct("beta",0.62,"p",0.04,"gamma",1/2.6,"N",62000000);
ICs8 = struct("S",Classes7.S(end),"Ir",Classes7.Ir(end),"In",Classes7.In(end),"R",Classes7.R(end));
mintime8 = 245;
maxtime8 = 500;
[Classes8] = ODE_SIR3_model(para8,ICs8,mintime8,maxtime8);

R08 = para8.beta/para8.gamma;
HIT8perc = 100*(1-(1/R08));

tim678 = [Classes6.t;Classes7.t;Classes8.t];
R678 = [Classes6.R;Classes7.R;Classes8.R];
%end of with initial immunity

%without initial immunity
para6w = struct("beta",0.62,"p",0.04,"gamma",1/2.6,"N",62000000);
ICs6w = struct("S",para6w.N-1,"Ir",0,"In",1,"R",0);
mintime6 = 119;
maxtime6 = 210;
[Classes6w] = ODE_SIR3_model(para6w,ICs6w,mintime6,maxtime6);

R06w = para6w.beta/para6.gamma;
HIT6percw = 100*(1-(1/R06));

para7w = struct("beta",0.62*0.6,"p",0.04,"gamma",1/2.6,"N",62000000);
ICs7w = struct("S",Classes6w.S(end),"Ir",Classes6w.Ir(end),"In",Classes6w.In(end),"R",Classes6w.R(end));
mintime7 = 210;
maxtime7 = 245;
[Classes7w] = ODE_SIR3_model(para7w,ICs7w,mintime7,maxtime7);

R07w = para7w.beta/para7w.gamma;
HIT7percw = 100*(1-(1/R07));

para8w = struct("beta",0.62,"p",0.04,"gamma",1/2.6,"N",62000000);
ICs8w = struct("S",Classes7w.S(end),"Ir",Classes7w.Ir(end),"In",Classes7w.In(end),"R",Classes7w.R(end));
mintime8 = 245;
maxtime8 = 500;
[Classes8w] = ODE_SIR3_model(para8w,ICs8w,mintime8,maxtime8);

R08w = para8w.beta/para8w.gamma;
HIT8percw = 100*(1-(1/R08));

tim678w = [Classes6w.t;Classes7w.t;Classes8w.t];
R678w = [Classes6w.R;Classes7w.R;Classes8w.R];
%end of without initial immunity

%without school holiday but with initial immunity

para9 = struct("beta",0.62,"p",0.04,"gamma",1/2.6,"N",62000000);
ICs9 = struct("S",0.87*para9.N-1,"Ir",0,"In",1,"R",0);
mintime9 = 119;
maxtime9 = 500;
[Classes9] = ODE_SIR3_model(para9,ICs9,mintime9,maxtime9);

R09 = para6.beta/para6.gamma;
HIT9perc = 100*(1-(1/R06));

%end of without school holiday

startt = 100*(Classes7.R(1)+0.13*para7.N)/para7.N;
endd = 100*(Classes8.R(1)+0.13*para7.N)/para7.N;

%plot hit and proportions
figure(16);
plot(tim678-119,100*(R678+0.13*para6.N)/para6.N);
hold on
plot(tim678w-119,100*(R678w)/para6.N);
hold on
plot(Classes9.t-119,100*(Classes9.R+0.13*para9.N)/para9.N);
hold on
yline(HIT6perc,"--");
hold on
plot(210-119,startt,"x","color",[0 0 255]./255);
hold on
plot(245-119,endd,"x","color",[255 0 0]./255);
hold off
xlabel("Time / days");
ylabel(["Percentage of N who are immune" ; "(R $+$ initial immunity $\%$"]);
legend("$13\%$ initial immunity with holiday","$0\%$ initial immunity with holiday","$13\%$ initial immunity without holiday","HIT $\approx 38 \%$","Start of holiday", "End of holiday",'Interpreter','latex');




%Q2 & 1 have same situation but with prior immunity


%paras and ics
para6 = struct("beta",0.62,"p",0.04,"gamma",1/2.6,"N",62000000);
ICs6 = struct("S",0.87*para6.N-1,"Ir",0,"In",1,"R",0);
maxtime6 = 250;
[Classes6] = ODE_SIR2_model(para6,ICs6,maxtime6);

para7 = struct("beta",0.62,"p",0.04,"gamma",1/2.6,"N",62000000);
ICs7 = struct("S",0.87*para7.N-1,"Ir",0,"In",1,"R",0);
maxtime7 = 250;
[Classes7] = ODE_SIR2_model(para7,ICs7,maxtime7);

%plot hit and proportions
figure(17);
plot(Classes6.t, 100*(Classes6.R+0.13*para6.N)/para6.N);
hold on
plot(Classes7.t, 100*Classes.R/para7.N)
hold on
yline(HIT6perc,"--");
hold off
xlabel("Time / days");
ylabel(["Percentage of N who are immune" ; "(R $+$ initial immunity $\%$"]);
legend("$13\%$ initial immunity","$0\%$ initial immunity","HIT $\approx 38 \%$",'Interpreter','latex');


