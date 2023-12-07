clear
clc

%Q3

%a
%import csv
tablg = readtable("SwineFluCaseData.csv");
dates = tablg{:,2};
cases = tablg{:,4};
cumcases = tablg{:,5};

%bar plot of weekly reported cases
figure(9);
bar(dates,cases);
xlabel("Week");
ylabel("Weekly Cases");
set(gca,"fontsize",16);
set(0,"defaultaxesfontsize",16);
set(0,"defaultlinelinewidth",2);
hold off

figure(10);
bar(dates,cumcases);
xlabel("Week");
ylabel("Weekly Cases Cumulative");
set(gca,"fontsize",16);
set(0,"defaultaxesfontsize",16);
set(0,"defaultlinelinewidth",2);


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
%new ics
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
plot(Classes3.t,Classes3.Ir+Classes3.In);
hold on
plot(Classes4.t,Classes4.Ir+Classes4.In);
hold on
plot(Classes5.t,Classes5.Ir+Classes5.In);
hold off
xlabel("Time / days");
ylabel("Number of infections, I(t)=Ir(t)+In(t)");
legend("Pre-holiday, \beta=0.62","During holiday, \beta=0.372","Post-holiday, \beta=0.62");

%c
Re3 = para3.beta.*Classes3.S./(para3.gamma*para3.N);
Re4 = para4.beta.*Classes4.S./(para4.gamma*para4.N);
Re5 = para5.beta.*Classes5.S./(para5.gamma*para5.N);

%plot
figure(12);
plot(Classes3.t,Re3);
hold on
plot(Classes4.t,Re4);
hold on
plot(Classes5.t,Re5);
hold on
xline(Classes3.t(end),":")
hold on
xline(Classes4.t(end),":")
hold off
xlabel("Time / days");
ylabel("Number of infections, I(t)=Ir(t)+In(t)");
legend("Pre-holiday, \beta=0.62","During holiday, \beta=0.372","Post-holiday, \beta=0.62");

%d

%final size by plot
fin_size = round(Classes5.R(end));

%Duration of outbreak, taken as when In+Ir<1
LastInft = find(Classes5.Ir+Classes5.In>=1,1,"last");
Duration = ceil(Classes5.t(LastInft));



%e 1
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
bar(dates,cases);
hold on
plot([17:70],wi);
hold off
xlabel("Time / weeks");
ylabel("Weekly reported incidence");
legend("Actual data", "Model Weekly incidence");



%e 2 by days scaled down
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

wit = nonzeros(wit);

%plot
figure(14);
bar(dates,cases);
hold on
plot(17:17+size(wit)-1,wit);
hold off
xlabel("Time / weeks");
ylabel("Weekly reported incidence vs data collected");
legend("Actual data", "Model Weekly incidence");


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
bar(dates,cases);
hold on
plot([17:70],ditot);
hold off
xlabel("Time / weeks");
ylabel("Weekly reported incidence vs data collected");
legend("Actual data", "Model Weekly incidence");



