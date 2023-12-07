E=0;
I=1;
R=0; 
D=0;
figure(100)
clf
S=10000-1;
[t,x1]=ode45('SEIRDfunction',[0 400],[S E I R D]);%ode45函数，用来计算微分方程
Read=[t,x1];
hold on %画图
line1=plot(t,x1(:,2),'color',[255/255 165/255 0/255],'linewidth',2);
line2=plot(t,x1(:,3),'color',[238/255 58/255 140/255],'linewidth',2);xmax=max(x1(:,3));yline(xmax,'-.b');
% id=find(x1(:,3)==xmax);ymax=t(id);xline(ymax,'-.b');annotation('textarrow',[0.5 0.123+ymax/100*0.81],[0.5 0.11+xmax/40000*0.81],'String','感染者人数达到极大值')
line3=plot(t,x1(:,5),'color',[171/255 130/255 255/255],'linewidth',2);
line4=plot(t,x1(:,4),'color',[135/255 206/255 255/255],'linewidth',2),legend([line1 line2 line3 line4],'Exposed','Infected','Death','Recovery');


% axis([0 100 0 80000]);%设置图表坐标
xlabel('Date'),ylabel('Person'),title('SEIRD Model','fontsize',13);%设置坐标轴和题目
hold off


figure(101);
plot(t,x1(:,3),"--k");