clear;
clc;
%%数据转换坐标系和消除零漂%%
n1=1;n2=12000;
A1=xlsread('ID3-1.xlsx',4);
%%测力仪x正方向和模型的下x正方向不一致
A1(:,3)=-A1(:,3);
% A1(:,2)=A1(:,2)-(sum(A1(1:1000,2))+sum(A1(m-999:m,2)))/2000;
% A1(:,3)=A1(:,3)-(sum(A1(1:1000,3))+sum(A1(m-999:m,3)))/2000;
% A1(:,4)=A1(:,4)-(sum(A1(1:1000,4))+sum(A1(m-999:m,4)))/2000;
figure(1)
subplot(3,1,1);
plot(A1(n1:n2,1),A1(n1:n2,2),'r');
grid on;
subplot(3,1,2);
plot(A1(n1:n2,1),A1(n1:n2,3),'b');
grid on;
subplot(3,1,3);
plot(A1(n1:n2,1),A1(n1:n2,4),'k');
grid on;
% subplot(3,1,1);
% plot(A1(:,1),A1(:,2),'r');
% grid on;
% subplot(3,1,2);
% plot(A1(:,1),A1(:,3),'b');
% grid on;
% subplot(3,1,3);
% plot(A1(:,1),A1(:,4),'k');
% grid on;

A2=xlsread('ID3-2.xlsx',4);
A2(:,3)=-A2(:,3);
% [m,n]=size(A2);
% A2(:,2)=A2(:,2)-(sum(A2(1:1000,2))+sum(A2(m-999:m,2)))/2000;
% A2(:,3)=A2(:,3)-(sum(A2(1:1000,3))+sum(A2(m-999:m,3)))/2000;
% A2(:,4)=A2(:,4)-(sum(A2(1:100,4))+sum(A2(m-99:m,4)))/200;
figure(2)
subplot(3,1,1);
plot(A2(n1:n2,1),A2(n1:n2,2),'r');
grid on;
subplot(3,1,2);
plot(A2(n1:n2,1),A2(n1:n2,3),'b');
grid on;
subplot(3,1,3);
plot(A2(n1:n2,1),A2(n1:n2,4),'k');
grid on;
% subplot(3,1,1);
% plot(A2(:,1),A2(:,2),'r');
% grid on;
% subplot(3,1,2);
% plot(A2(:,1),A2(:,3),'b');
% grid on;
% subplot(3,1,3);
% plot(A2(:,1),A2(:,4),'k');
% grid on;

A3=xlsread('ID3-3.xlsx',4);
A3(:,3)=-A3(:,3);
% A3(:,2)=A3(:,2)-(sum(A3(1:1000,2))+sum(A3(m-999:m,2)))/2000;
% A3(:,3)=A3(:,3)-(sum(A3(1:1000,3))+sum(A3(m-999:m,3)))/2000;
% A3(:,4)=A3(:,4)-(sum(A3(1:100,4))+sum(A3(m-99:m,4)))/200;
figure(3)
subplot(3,1,1);
plot(A3(n1:n2,1),A3(n1:n2,2),'r');
grid on;
subplot(3,1,2);
plot(A3(n1:n2,1),A3(n1:n2,3),'b');
grid on;
subplot(3,1,3);
plot(A3(n1:n2,1),A3(n1:n2,4),'k');
grid on;
% subplot(3,1,1);
% plot(A3(:,1),A3(:,2),'r');
% grid on;
% subplot(3,1,2);
% plot(A3(:,1),A3(:,3),'b');
% grid on;
% subplot(3,1,3);
% plot(A3(:,1),A3(:,4),'k');
% grid on;

A4=xlsread('ID3-4.xlsx',4);
A4(:,3)=-A4(:,3);
% A4(:,2)=A4(:,2)-(sum(A4(1:1000,2))+sum(A4(m-999:m,2)))/2000-5;
% A4(:,3)=A4(:,3)-(sum(A4(1:1000,3))+sum(A4(m-999:m,3)))/2000;
% A4(:,4)=A4(:,4)-(sum(A4(1:100,4))+sum(A4(m-99:m,4)))/200;
figure(4)
subplot(3,1,1);
plot(A4(n1:n2,1),A4(n1:n2,2),'r');
grid on;
subplot(3,1,2);
plot(A4(n1:n2,1),A4(n1:n2,3),'b');
grid on;
subplot(3,1,3);
plot(A4(n1:n2,1),A4(n1:n2,4),'k');
grid on;
% subplot(3,1,1);
% plot(A4(:,1),A4(:,2),'r');
% grid on;
% subplot(3,1,2);
% plot(A4(:,1),A4(:,3),'b');
% grid on;
% subplot(3,1,3);
% plot(A4(:,1),A4(:,4),'k');
% grid on;

A5=xlsread('ID3-5.xlsx',4);
A5(:,3)=-A5(:,3);
% A5(:,2)=A5(:,2)-(sum(A5(1:1000,2))+sum(A5(m-999:m,2)))/2000-1;
% A5(:,3)=A5(:,3)-(sum(A5(1:1000,3))+sum(A5(m-999:m,3)))/2000;
% A5(:,4)=A5(:,4)-(sum(A5(1:100,4))+sum(A5(m-99:m,4)))/200;
figure(5)
subplot(3,1,1);
plot(A5(n1:n2,1),A5(n1:n2,2),'r');
grid on;
subplot(3,1,2);
plot(A5(n1:n2,1),A5(n1:n2,3),'b');
grid on;
subplot(3,1,3);
plot(A5(n1:n2,1),A5(n1:n2,4),'k');
grid on;
% subplot(3,1,1);
% plot(A5(:,1),A5(:,2),'r');
% grid on;
% subplot(3,1,2);
% plot(A5(:,1),A5(:,3),'b');
% grid on;
% subplot(3,1,3);
% plot(A5(:,1),A5(:,4),'k');
% grid on;

Av=zeros(5,4);
Av(1,2)=mean(A1(n1:n2,2));Av(1,3)=mean(A1(n1:n2,3));Av(1,4)=mean(A1(n1:n2,4));
Av(2,2)=mean(A2(n1:n2,2));Av(2,3)=mean(A2(n1:n2,3));Av(2,4)=mean(A2(n1:n2,4));
Av(3,2)=mean(A3(n1:n2,2));Av(3,3)=mean(A3(n1:n2,3));Av(3,4)=mean(A3(n1:n2,4));
Av(4,2)=mean(A4(n1:n2,2));Av(4,3)=mean(A4(n1:n2,3));Av(4,4)=mean(A4(n1:n2,4));
Av(5,2)=mean(A5(n1:n2,2));Av(5,3)=mean(A5(n1:n2,3));Av(5,4)=mean(A5(n1:n2,4));
Av(1,1)=0.03;Av(2,1)=0.04;Av(3,1)=0.06;Av(4,1)=0.08;Av(5,1)=0.10;

yx=polyfit(Av(:,1),Av(:,2),1);
y1=polyval(yx,0:0.01:0.1);
yy=polyfit(Av(:,1),Av(:,3),1);
y2=polyval(yy,0:0.01:0.1);
yz=polyfit(Av(:,1),Av(:,4),1);
y3=polyval(yz,0:0.01:0.1);
figure(6)
plot(Av(:,1),Av(:,2),'ro');
hold on;
plot(Av(:,1),Av(:,3),'bo');
hold on;
plot(Av(:,1),Av(:,4),'ko');
hold on;
plot(0:0.01:0.1,y1,'r-');
hold on;
plot(0:0.01:0.1,y2,'b-');
hold on;
plot(0:0.01:0.1,y3,'k-');
hold on;
title('Average forces in different experiments');
legend('Mean force in x','Mean force in y','Mean force in z');
ap=0.6;%%切深0.6mm
N=2;%%齿数为1
Ktc=4*(yy(1))/(N*ap);Kte=pi*(yy(2))/(N*ap);
Krc=-4*(yx(1))/(N*ap);Kre=-pi*(yx(2))/(N*ap);
Kac=pi*(yz(1))/(N*ap);Kae=2*(yz(2))/(N*ap);