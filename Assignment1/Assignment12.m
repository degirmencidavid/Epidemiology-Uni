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
legend("S","I","R");
set(gca,"fontsize",16);
set(0,"defaultaxesfontsize",16);
set(0,"defaultlinelinewidth",2);


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

%daily incidence
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
set(gca,"fontsize",16);
set(0,"defaultaxesfontsize",16);
set(0,"defaultlinelinewidth",2);

%plot cumulative number of reported infections
figure(7);
plot(Classes2.t,round(trc))
xlabel("Time / days");
ylabel("Total Reported Cases");
set(gca,"fontsize",16);
set(0,"defaultaxesfontsize",16);
set(0,"defaultlinelinewidth",2);

%check if cumulative reported cases makes sense (they do :), very nice)
%disp(round(trc(end)*25)-round(Classes2.R(end)))

















% 
% 
% %cumulative
% FOI2 = para2.beta.*Classes2.S.*(Classes2.Ir+Classes2.In)./para2.N;
% NewInfections = round(trapz(Classes2.t,para2.p.*FOI2));
% Final_size_reported_I = NewInfections + para2.p.*(ICs2.Ir+ICs2.In);
% Final_size_I = NewInfections + ICs2.Ir+ICs2.In;
% trc = cumtrapz(Classes2.t,para2.p.*FOI2);
% 
% %plot cumulative
% figure(6);
% plot(Classes2.t, trc)
% xlabel("Time / days");
% ylabel("Total Reported Cases");
% set(gca,"fontsize",16);
% set(0,"defaultaxesfontsize",16);
% set(0,"defaultlinelinewidth",2);
% 
% %daily incidence
% dit(1) = Classes2.Ir(1);
% for i=2:floor(size(trc))
%     dit(i) = trc(i) - trc(i-1);
% end
% 
% %dit = trapz(Classes2.t,para2.p./FOI2);
% 
% %plot daily incidence probably a mistake somewhere
% figure(7);
% plot(Classes2.t,dit);
% xlabel("Time / days");
% ylabel("Daily Reported Incidence");
% set(gca,"fontsize",16);
% set(0,"defaultaxesfontsize",16);
% set(0,"defaultlinelinewidth",2);
% 
% 
% 
% %fig 8 prolly wrong
% %daily
% %Incidence is the number of new infections between two time points
%  for it=1:floor(size(FOI2)-1)
%     di(it) = trapz(Classes2.t([(it-1):it]+1), para2.p.*FOI2([(it-1):it]+1));
%  end
%  
% %  for it=1:floor(size(FOI2)-1)
% %     di(it) = trapz(Classes.t([(it-1):it]+1), FOI2([(it-1):it]+1));
% %  end
%  
%  
% %plot daily incidence probably a mistake somewhere
% figure(8);
% plot([1:size(FOI2)-1], para2.p.*di);
% xlabel("Time / days");
% ylabel("Daily Reported Incidence");
% set(gca,"fontsize",16);
% set(0,"defaultaxesfontsize",16);
% set(0,"defaultlinelinewidth",2);
% 
