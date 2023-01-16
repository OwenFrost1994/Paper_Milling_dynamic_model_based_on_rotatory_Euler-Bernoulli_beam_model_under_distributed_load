clear;
clc;
close all;
%%专门计算梁的模态振型
%%首先生成梁模型的基本参数，计算矩阵%%
%%刀具等效为定截面两段静态梁模型%%
DL=0.02;%%梁单元的长度,mm
E=577500000000;%%杨氏模量N/m^2
p=14500;%%密度kg/(m^3)
L0=60;%%梁总悬长，mm，0号节点不用管，值已知
L0=L0/1000;%%梁总悬长，m
L1=30;%%第一段截面的长度，mm
L1=L1/1000;%%第一段截面的长度，m
L2=L0-L1;%%第二段截面的长度，mm
D1=6;%%第一段直径，mm
D1=D1/1000;%%第一段直径，m
D2=4.8;%%第二段直径，mm
D2=D2/1000;%%第二段直径，m
I1=pi*(D1)^4/64;%%第一段惯性矩,m^4
I2=pi*(D2)^4/64;%%第二段惯性矩,m^4
DL=DL/1000;%%梁单元的长度,m
NB=ceil(L0/DL);%%梁单元的个数

%%惯性矩离散%%
for i=1:1:NB
    l(i)=i*DL-0.5*DL;
    if l(i)<=L1
        I(i)=I1;
    else
        I(i)=I2;
    end
end

%%截面积离散%%
%%这里的i是从1到N一共N个节点，没有在原有节点的基础上加虚拟节点
for i=1:1:NB
    l(i)=i*DL-0.5*DL;%%节点的位置处于每个片层的中点
    if l(i)<=L1%%前半段截面积
        S(i)=pi*D1^2/4;
    else
        S(i)=pi*D2^2/4;%%后半段截面积
    end
end

%%计算I'',I',I，存入ID中，每一个节点的这个三个值都是行向量
for i=1:1:NB
    if i==1
        I2=(I(i)-2*I(i+1)+I(i+2))/(DL^2);
        I1=(-3*I(i)+4*I(i+1)-I(i+2))/(2*DL);
    else if i==NB
                    I2=(I(i-2)-2*I(i-1)+I(i))/(DL^2);
                    I1=(I(i-2)-4*I(i-1)+3*I(i))/(2*DL);
        else
            I2=(I(i-1)-2*I(i)+I(i+1))/(DL^2);
            I1=(I(i+1)-I(i-1))/(2*DL);
        end
    end
    ID(i,:)=[I2,2*I1,I(i)];
end

%%以上为准备工作%%
%%下面根据每个时刻生成系数矩阵，常数矩阵，和矩阵求解%%
%%这里前后虚拟节点已经在公式当中考虑过了%%
%%产生系数矩阵%%
A1=zeros(NB,NB);
A2=zeros(NB,NB);
A3=zeros(NB,NB);
%%二阶导的系数矩阵%%
for i=1:1:NB
    if i==1
        A1(i,1:5)=[-31,16,-1,0,0];
    else if i==2
            A1(i,1:5)=[16,-30,16,-1,0];
        else if i==NB-1
                A1(i,NB-4:NB)=[0,-1,111/7,-201/7,97/7];
            else if i==NB
                    A1(i,NB-4:NB)=[0,0,0,0,0];
                else
                    A1(i,i-2:i+2)=[-1,16,-30,16,-1];
                end
            end
        end
    end
end
%%三阶导的系数矩阵%%
for i=1:1:NB
    if i==1
        A2(i,1:5)=[-1,-2,1,0,0];
    else if i==2
            A2(i,1:5)=[2,0,-2,1,0];
        else if i==NB-1
                A2(i,NB-4:NB)=[0,-1,15/7,-9/7,1/7];
            else if i==NB
                    A2(i,NB-4:NB)=[0,0,0,0,0];
                else
                    A2(i,i-2:i+2)=[-1,2,0,2,-1];
                end
            end
        end
    end
end
%%四阶导的系数矩阵%%
for i=1:1:NB
    if i==1
        A3(i,1:5)=[7,-4,1,0,0];
    else if i==2
            A3(i,1:5)=[-4,6,-4,1,0];
        else if i==NB-1
                A3(i,NB-4:NB)=[0,1,-27/7,33/7,-13/7];
            else if i==NB
                    A3(i,NB-4:NB)=[0,0,12/7,-24/7,12/7];
                else
                    A3(i,i-2:i+2)=[1,-4,6,-4,1];
                end
            end
        end
    end
end
%%获得刚度矩阵%%
for i=1:1:NB
    KS(i,:)=E*DL*(ID(i,1)*A1(i,:)/(12*DL^2)+ID(i,2)*A2(i,:)/(2*DL^3)+ID(i,3)*A3(i,:)/(DL^4));
end
%%获得质量矩阵%%
for i=1:1:NB
    M(i,i)=p*DL*S(i);
end

DL=DL*1000;
L0=L0*1000;
%%mesh函数绘制扁平梁
[U,Wn]=eigs(inv(M)*KS,6,'SM');
Wn=sqrt(Wn)/(2*pi);

figure(1)
XA=zeros(NB,1);
X=[XA-0.1,XA,XA+0.1];
Y=[DL:DL:L0;DL:DL:L0;DL:DL:L0]';
Z=[U(:,1),U(:,1),U(:,1)];
mesh(X,Y,Z)
% plot3(XA,DL:DL:L0,U(:,1),'k','linewidth',6.0)
% patch([XA' nan],[DL:DL:L0 nan],[U(:,1)' nan],'edgecolor','interp')
hold on;
XA=1*ones(NB,1);
X=[XA-0.1,XA,XA+0.1];
Z=[U(:,2),U(:,2),U(:,2)];
mesh(X,Y,Z)
% plot3(XA,DL:DL:L0,U(:,2),'r','linewidth',6.0)%r
hold on;
XA=2*ones(NB,1);
X=[XA-0.1,XA,XA+0.1];
Z=[U(:,3),U(:,3),U(:,3)];
mesh(X,Y,Z)
% plot3(XA,DL:DL:L0,U(:,3),'b','linewidth',6.0)%b
hold on;
XA=3*ones(NB,1);
X=[XA-0.1,XA,XA+0.1];
Z=[U(:,4),U(:,4),U(:,4)];
mesh(X,Y,Z)
% plot3(XA,DL:DL:L0,U(:,4),'y','linewidth',6.0)%y
hold on;
XA=4*ones(NB,1);
X=[XA-0.1,XA,XA+0.1];
Z=[U(:,5),U(:,5),U(:,5)];
mesh(X,Y,Z)
% plot3(XA,DL:DL:L0,U(:,5),'g','linewidth',6.0)%g
hold on;
XA=5*ones(NB,1);
X=[XA-0.1,XA,XA+0.1];
Z=[U(:,6),U(:,6),U(:,6)];
mesh(X,Y,Z)
% plot3(XA,DL:DL:L0,U(:,6),'m','linewidth',6.0)%m
grid on;
hold on;
ylim([0 L0])
zlim([-0.1 0.1])
xlabel('ω order')
ylabel('L(mm)')
zlabel('Displacement(mm)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
set(get(gca,'ZLabel'),'Fontsize',20)
legend(strcat('ω_0=',num2str(Wn(1,1)),'Hz'),strcat('ω_1=',num2str(Wn(2,2)),'Hz'),strcat('ω_2=',num2str(Wn(3,3)),'Hz'),strcat('ω_3=',num2str(Wn(4,4)),'Hz'),strcat('ω_4=',num2str(Wn(5,5)),'Hz'),strcat('ω_5=',num2str(Wn(6,6)),'Hz'))
view(27.5, 60);

SS=1000;%%主轴转速rpm
w=2*pi*SS/60;%%刀具角速度
[i,j]=size(KS);
KW=KS-w^2*p*diag(S);

[Uw,Wnw]=eigs(inv(M)*KW,6,'SM');
Wnw=sqrt(Wnw)/(2*pi);

figure(2)
XA=zeros(NB,1);
X=[XA-0.1,XA,XA+0.1];
Y=[DL:DL:L0;DL:DL:L0;DL:DL:L0]';
Z=[Uw(:,1),Uw(:,1),Uw(:,1)];
mesh(X,Y,Z)
hold on;
XA=1*ones(NB,1);
X=[XA-0.1,XA,XA+0.1];
Z=[Uw(:,2),Uw(:,2),Uw(:,2)];
mesh(X,Y,Z)
hold on;
XA=2*ones(NB,1);
X=[XA-0.1,XA,XA+0.1];
Z=[Uw(:,3),Uw(:,3),Uw(:,3)];
mesh(X,Y,Z)
hold on;
XA=3*ones(NB,1);
X=[XA-0.1,XA,XA+0.1];
Z=[Uw(:,4),Uw(:,4),Uw(:,4)];
mesh(X,Y,Z)
hold on;
XA=4*ones(NB,1);
X=[XA-0.1,XA,XA+0.1];
Z=[Uw(:,5),Uw(:,5),Uw(:,5)];
mesh(X,Y,Z)
hold on;
XA=5*ones(NB,1);
X=[XA-0.1,XA,XA+0.1];
Z=[Uw(:,6),Uw(:,6),Uw(:,6)];
mesh(X,Y,Z)
grid on;
hold on;
ylim([0 L0])
zlim([-0.1 0.1])
xlabel('ω order')
ylabel('L(mm)')
zlabel('Displacement(mm)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
set(get(gca,'ZLabel'),'Fontsize',20)
legend(strcat('ω_0=',num2str(Wnw(1,1)),'Hz'),strcat('ω_1=',num2str(Wnw(2,2)),'Hz'),strcat('ω_2=',num2str(Wnw(3,3)),'Hz'),strcat('ω_3=',num2str(Wnw(4,4)),'Hz'),strcat('ω_4=',num2str(Wnw(5,5)),'Hz'),strcat('ω_5=',num2str(Wnw(6,6)),'Hz'))
view(27.5, 60);

figure(3)
subplot(3,2,1)
plot(DL:DL:L0,Uw(:,1),'r-','Markersize',7,'Markerface','white','linewidth',6.0)
hold on;
plot(DL:DL:L0,U(:,1),'b-.','Markersize',7,'Markerface','white','linewidth',6.0)
hold on;
grid on;
xlabel('L(mm)')
ylabel('Displacement(mm)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
legend(strcat('ω_r0=',num2str(Wnw(1,1)),'Hz'),strcat('ω_s0=',num2str(Wn(1,1)),'Hz'))

subplot(3,2,2)
plot(DL:DL:L0,Uw(:,2),'r-','Markersize',7,'Markerface','white','linewidth',6.0)
hold on;
plot(DL:DL:L0,U(:,2),'b-.','Markersize',7,'Markerface','white','linewidth',6.0)
hold on;
grid on;
xlabel('L(mm)')
ylabel('Displacement(mm)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
legend(strcat('ω_r1=',num2str(Wnw(2,2)),'Hz'),strcat('ω_s1=',num2str(Wn(2,2)),'Hz'))

subplot(3,2,3)
plot(DL:DL:L0,Uw(:,3),'r-','Markersize',7,'Markerface','white','linewidth',6.0)
hold on;
plot(DL:DL:L0,U(:,3),'b-.','Markersize',7,'Markerface','white','linewidth',6.0)
hold on;
grid on;
xlabel('L(mm)')
ylabel('Displacement(mm)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
legend(strcat('ω_r2=',num2str(Wnw(3,3)),'Hz'),strcat('ω_s2=',num2str(Wn(3,3)),'Hz'))

subplot(3,2,4)
plot(DL:DL:L0,Uw(:,4),'r-','Markersize',7,'Markerface','white','linewidth',6.0)
hold on;
plot(DL:DL:L0,U(:,4),'b-.','Markersize',7,'Markerface','white','linewidth',6.0)
hold on;
grid on;
xlabel('L(mm)')
ylabel('Displacement(mm)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
legend(strcat('ω_r3=',num2str(Wnw(4,4)),'Hz'),strcat('ω_s3=',num2str(Wn(4,4)),'Hz'))

subplot(3,2,5)
plot(DL:DL:L0,Uw(:,5),'r-','Markersize',7,'Markerface','white','linewidth',6.0)
hold on;
plot(DL:DL:L0,U(:,5),'b-.','Markersize',7,'Markerface','white','linewidth',6.0)
hold on;
grid on;
xlabel('L(mm)')
ylabel('Displacement(mm)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
legend(strcat('ω_r4=',num2str(Wnw(5,5)),'Hz'),strcat('ω_s4=',num2str(Wn(5,5)),'Hz'))

subplot(3,2,6)
plot(DL:DL:L0,Uw(:,6),'r-','Markersize',7,'Markerface','white','linewidth',6.0)
hold on;
plot(DL:DL:L0,U(:,6),'b-.','Markersize',7,'Markerface','white','linewidth',6.0)
hold on;
grid on;
xlabel('L(mm)')
ylabel('Displacement(mm)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
legend(strcat('ω_r5=',num2str(Wnw(6,6)),'Hz'),strcat('ω_s5=',num2str(Wn(6,6)),'Hz'))