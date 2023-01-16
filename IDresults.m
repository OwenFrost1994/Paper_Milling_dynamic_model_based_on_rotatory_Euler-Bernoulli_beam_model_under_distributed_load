clear;
clc;
Ktc=xlsread('coefficientswithAE.xlsx',1);
Kte=xlsread('coefficientswithAE.xlsx',2);
Krc=xlsread('coefficientswithAE.xlsx',3);
Kre=xlsread('coefficientswithAE.xlsx',4);
Ktcave=xlsread('coefficientswithAE.xlsx',5);
Kteave=xlsread('coefficientswithAE.xlsx',6);
Krcave=xlsread('coefficientswithAE.xlsx',7);
Kreave=xlsread('coefficientswithAE.xlsx',8);

pKtc=polyfit(Ktcave(1,:),Ktcave(2,:),3);
yKtc=polyval(pKtc,0:0.1:6);
pKte=polyfit(Kteave(1,:),Kteave(2,:),3);
yKte=polyval(pKte,0:0.1:6);
pKrc=polyfit(Krcave(1,:),Krcave(2,:),3);
yKrc=polyval(pKrc,0:0.1:6);
pKre=polyfit(Kreave(1,:),Kreave(2,:),3);
yKre=polyval(pKre,0:0.1:6);

figure(1)

subplot(2,2,1)
plot(0:0.1:6,yKtc,'r-.','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
plot(Ktc(1,:),Ktc(2,:),'b+','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
grid on;
xlabel('ae(mm)')
ylabel('Ktc(N/mm^2)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('Cuttig Force Coefficients with Radial Depth of Cut');
subplot(2,2,2)
plot(0:0.1:6,yKrc,'r-.','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
plot(Krc(1,:),Krc(2,:),'b+','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
grid on;
xlabel('ae(mm)')
ylabel('Krc(N/mm^2)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
subplot(2,2,3)
plot(0:0.1:6,yKte,'r-.','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
plot(Kte(1,:),Kte(2,:),'b+','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
grid on;
xlabel('ae(mm)')
ylabel('Kte(N/mm)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
subplot(2,2,4)
plot(0:0.1:6,yKre,'r-.','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
plot(Kre(1,:),Kre(2,:),'b+','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
grid on;
xlabel('ae(mm)')
ylabel('Kre(N/mm)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)

