clear;
clc;
close all;
%%梁动力学系统在切削力下的仿真，刀具等效为两段静态梁模型%%首先生成梁模型的基本参数，计算矩阵%%
DL=0.2;%%梁单元的长度,mm
E=719000000000;%%杨氏模量N/m^2
p=14400;%%密度kg/(m^3)
L0=60;%%梁总悬长，mm
L1=30;%%第一段截面的长度，mm
L2=L0-L1;%%第二段截面的长度，mm
D1=6;%%第一段直径，mm
D2=4.8;%%第二段直径，mm
I1=pi*(D1/1000)^4/64;%%第一段惯性矩,m^4
S1=pi*(D1/1000)^2/4;%%第一段截面积,m^2
I2=pi*(D2/1000)^4/64;%%第二段惯性矩,m^4
S2=pi*(D2/1000)^2/4;%%第二段截面积,m^2
NB=ceil(L0/DL);%%梁单元的个数

%%惯性矩和截面积离散%%
for i=1:1:NB
    l(i)=i*DL-0.5*DL;%%计算层片的中点位置
    if l(i)<=L1
        I(i)=I1;
        S(1,i)=S1;
    else
        I(i)=I2;
        S(1,i)=S2;
    end
end

%%计算I'',I',I，存入ID中，每一个节点的这个三个值都是行向量
for i=1:1:NB
    if i==1
        I2=(I(i)-2*I(i+1)+I(i+2))/((0.001*DL)^2);
        I1=(-3*I(i)+4*I(i+1)-I(i+2))/(2*(0.001*DL));
    else if i==NB
                    I2=(I(i-2)-2*I(i-1)+I(i))/((0.001*DL)^2);
                    I1=(I(i-2)-4*I(i-1)+3*I(i))/(2*(0.001*DL));
        else
            I2=(I(i-1)-2*I(i)+I(i+1))/((0.001*DL)^2);
            I1=(I(i+1)-I(i-1))/(2*(0.001*DL));
        end
    end
    ID(i,:)=[I2,2*I1,I(i)];
end

%%以上为准备工作%%
%%下面根据每个时刻生成系数矩阵，常数矩阵，最终获得静态刚度矩阵和质量矩阵%%
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
    KS(i,:)=E*DL*(ID(i,1)*A1(i,:)/(12*(DL/1000)^2)+ID(i,2)*A2(i,:)/(2*(DL/1000)^3)+ID(i,3)*A3(i,:)/((DL/1000)^4))/1000;
end
%%获得质量矩阵%%
for i=1:1:NB
    M(i,i)=p*DL*S(i)/1000;
end

%%加入刀具参数，不考虑偏心和倾斜

%%刀具几何参数%%
D=D1;%%刀具半径
N=2;%%刀具齿数
B=pi/6;%%刀具螺旋角
Cp=2*pi/N;%%齿间角
%%微元长度即为上面梁单元的微元长度

%%材料参数%%材料AL6061-T6511%%
Ktc=613.92;%%切向剪切力系数
Kte=12.785;%%切向刃口力系数
Krc=330.75;%%径向剪切力系数
Kre=12.902;%%径向刃口力系数
Kac=78.03;%%轴向剪切力系数
Kae=1.3672;%%轴向刃口力系数

%%加工参数%%
%%铣削方式:顺铣%%
Cm=1;%%铣削方式，顺铣为1，逆铣为0
SS=2000;%%主轴转速
f=400;%%进给速度
fs=10000;%%采样频率
ap=4;%%轴向切深（单位mm）
ae=1*D;%%径向切深（单位mm）
Cn=5;%%圈数circle number

%%基本参数计算%%
R=D/2;%%刀具半径
kb=(2*tan(B))/D;%%kβ计算
fe=f/(N*SS);%%feed every tooth
w=2*pi*SS/60;%%刀具角速度
T=2*pi/w;%%刀具周期
Ns=floor(60*fs/SS);%%一个周期内的采样点个数
if Cm==1%%顺铣
    Cst=pi-acos((R-ae)/R);%%切入角
    Cex=pi;%%切出角
else%%逆铣
    Cst=0;%%切入角
    Cex=acos((R-ae)/R);%%切出角
end
Cs=0;%%开始角度
Dt=T/Ns;%%时间步长
DC=Dt*w;%%角度增量
Ca=Cs;%%初始角度

%%计算广义动刚度矩阵和广义阻尼矩阵，建立广义刚度和广义质量矩阵
[i,j]=size(KS);
KW=-w^2*p*diag([S,S]);
CQ=[zeros(i,j),diag(-2*w*p*S*DL/1000);diag(2*w*p*S*DL/1000),zeros(i,j)];
KQ=[KS,zeros(i,j);zeros(i,j),KS]+KW;
MQ=[M,zeros(i,j);zeros(i,j),M];

%%各种存储单元%%
Dl=DL;
FX=zeros(NB,Ns*Cn);
FY=zeros(NB,Ns*Cn);
FZ=zeros(NB,Ns*Cn);%%记录每一个梁单元上每一时刻的受力，其实际上只有头部的梁单元受力
Fx=0;
Fy=0;
Fz=0;
F=zeros(3,Ns*Cn);%%存储三个方向的切削力
x=zeros(NB,Ns*Cn);
y=zeros(NB,Ns*Cn);
Dx=zeros(ap/Dl,Ns*Cn);%%刀具在x方向的位移变化量，fv的第一个参数,这里动态切屑厚度只考虑切削段的部分，编号方式和刃编号方式一致
Dy=zeros(ap/Dl,Ns*Cn);%%刀具在y方向的位移变化量，fv的第二个参数
apx=zeros(1,Cn*N);
apy=zeros(1,Cn*N);

%%微分方程离散化求解
K1=MQ/(Dt^2)+3*CQ/(2*Dt)+KQ;
K2=2*MQ/(Dt^2)+2*CQ/Dt;
K3=-MQ/(Dt^2)-CQ/2*Dt;
K1=inv(K1);
%%计算整个模拟域的大循环
%%计算切削力同步用有限差分法计算梁的振动位移和速度%%
for i=1:1:Ns*Cn;
    Ca=Ca+DC;%%微元角度叠加计算刀具的转动角
    if Ca>=2*pi%%考虑刀具多个旋转周期，累加的刀具角度超过一周就减去一个2π
        Ca=Ca-2*pi;
    else
    end
    %%考虑动态切屑厚度%%
    for m=1:1:ap/Dl
        if i<=(Ns/N)%%如果是第一个刀齿在切的时候，切屑上表面没有前一个刀齿切出的波纹表面，此时相当于刀具无振动的情况
            Dx(m,i)=1000*(x(NB-m+1,i)-0);%%注意到动力学方程求解过程中单位为国际制单位，求解获得的位移单位为m，而每齿进给量的单位是mm
            Dy(m,i)=1000*(y(NB-m+1,i)-0);%%注意时间是每一列，不同刀刃编号是行
        else
            q=floor(i*N/Ns);
            if (i-q*Ns/N)==0
                q=q-1;%%如果正好是某个刀齿交替阶段，实际上需要减去一个刀齿
            else
            end
            for n=1:1:q%%循环取出前面的节点
                apx(n)=x(NB-m+1,i-n*Ns/N);%%注意沿刀轴方向的每个刀齿都要计算
                apy(n)=y(NB-m+1,i-n*Ns/N);%%注意沿刀轴方向的每个刀齿都要计算
            end
            Dx(m,i)=1000*(x(i)-min(apx(1:q)));
            Dy(m,i)=1000*(y(i)-min(apy(1:q)));
        end
    end

    %%最终将动态切屑厚度存入矩阵当中%%
    %%微元叠加计算一个刀刃上的切削力
    for m=1:1:ap/Dl
        Cd=Ca-m*kb*Dl+0.5*kb*Dl;
        %叠加计算多个刀齿的内循环
        for j=1:1:N
        C=Cd-(j-1)*Cp;%%考虑多齿存在的齿间角滞后
        if C<0
            C=C+2*pi;%%如果刀齿角度小于零则转成正的角度
        else
        end
        fa=fe*sin(C)++Dx(m,i)*sin(C)+Dy(m,i)*cos(C);%未变形切屑厚度
        if C<Cex&&C>Cst&&fa>=0
            fx=(-cos(C))*(Ktc*fa+Kte)*Dl+(-sin(C))*(Krc*fa+Kre)*Dl;
            fy=( sin(C))*(Ktc*fa+Kte)*Dl+(-cos(C))*(Krc*fa+Kre)*Dl;
            fz=(Kac*fa+Kae)*Dl;
        else
            fx=0;fy=0;fz=0;
        end
        Fx=Fx+fx;Fy=Fy+fy;Fz=Fz+fz;%%累加每一个刀片上每个齿所受的切削力
        end
        FX(NB-m+1,i)=Fx;FY(NB-m+1,i)=Fy;FZ(NB-m+1,i)=Fz;
        Fx=0;Fy=0;Fz=0;%%刀具受力累加变量归零
    end
    F(1,i)=sum(FX(:,i));F(2,i)=sum(FY(:,i));F(3,i)=sum(FZ(:,i));%%在矩阵中存储切削力
    %%xy方向的振动需要同时计算
    if  i==1 %%判断是否是初始点（第一个点）
        Q1=zeros(2*NB,1);
        Q2=zeros(2*NB,1);
    else if i==2%%不是初始点
            Q1=[x(:,i-1);y(:,i-1)];
            Q2=zeros(2*NB,1); 
        else
            Q1=[x(:,i-1);y(:,i-1)];
            Q2=[x(:,i-2);y(:,i-2)];
        end
    end
    Q=K1*([FX(:,i);FY(:,i)]+K2*Q1+K3*Q2);
    x(:,i)=Q(1:NB);
    y(:,i)=Q(NB+1:2*NB);
end

%%绘制力的图像，蓝色X向，红色Y向，绿色Z向
figure(1)
subplot(3,1,1)
plot(Dt:Dt:Cn*Ns*Dt,F(1,:),'b-','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
grid on;
xlabel('time(s)')
ylabel('Predicted F_x(N)')
set(gca, 'FontName','Times New Roman','FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('Forces imposed on cutter in three directions');
subplot(3,1,2)
plot(Dt:Dt:Cn*Ns*Dt,F(2,:),'r-','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
grid on;
xlabel('time(s)')
ylabel('Predicted F_y(N)')
set(gca, 'FontName','Times New Roman','FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
subplot(3,1,3)
plot(Dt:Dt:Cn*Ns*Dt,F(3,:),'k-','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
grid on;
xlabel('time(s)')
ylabel('Predicted F_z(N)')
set(gca, 'FontName','Times New Roman','FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)

figure(2)
plot(Dt:Dt:Cn*Ns*Dt,F(1,:),'b-','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
plot(Dt:Dt:Cn*Ns*Dt,F(2,:),'r-','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
plot(Dt:Dt:Cn*Ns*Dt,F(3,:),'k-','Markersize',7,'Markerface','white','linewidth',3.0);
grid on;
xlabel('time(s)')
ylabel('Force(N)')
set(gca, 'FontName','Times New Roman','FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('Summary of forces imposed on cutter');
legend('X predicted','Y predicted','Z predicted')

figure(3)
surf(Dt:Dt:Dt*Ns*Cn,DL:DL:L0,x)
shading FLAT;
colormap hot;
grid on;
hold on;
xlabel('time(s)')
ylabel('Length(m)')
zlabel('Displacement in x direction(m)')
set(gca, 'FontName','Times New Roman','FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
set(get(gca,'ZLabel'),'Fontsize',20)
title('Vibration displacement in x');
figure(4)
surf(Dt:Dt:Dt*Ns*Cn,DL:DL:L0,y)
shading FLAT;
colormap cool;
grid on;
hold on;
xlabel('time(s)')
ylabel('Length(m)')
zlabel('Displacement in y direction(m)')
set(gca, 'FontName','Times New Roman','FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
set(get(gca,'ZLabel'),'Fontsize',20)
title('Vibration displacement in y');
% 
% figure(5)
% surf(Dt:Dt:Dt*Ns*Cn,Dl:Dl:ap,Dx)
% grid on;
% hold on;
% xlabel('time(s)')
% ylabel('Length(m)')
% zlabel('Dynamic displacement difference in x direction(mm)')
% set(gca, 'FontName','Times New Roman','FontSize', 20)
% set(get(gca,'XLabel'),'Fontsize',20)
% set(get(gca,'YLabel'),'Fontsize',20)
% set(get(gca,'ZLabel'),'Fontsize',20)
% title('Dynamic displacement difference in x');
% 
% figure(6)
% surf(Dt:Dt:Dt*Ns*Cn,Dl:Dl:ap,Dy)
% grid on;
% hold on;
% xlabel('time(s)')
% ylabel('Length(m)')
% zlabel('Dynamic displacement difference in y direction(mm)')
% set(gca, 'FontName','Times New Roman','FontSize', 20)
% set(get(gca,'XLabel'),'Fontsize',20)
% set(get(gca,'YLabel'),'Fontsize',20)
% set(get(gca,'ZLabel'),'Fontsize',20)
% title('Dynamic displacement difference in y');
% 
figure(7)
subplot(3,1,1)
plot(Dt:Dt:Dt*Ns*Cn,x(NB,:),'r+-','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
plot(Dt:Dt:Dt*Ns*Cn,y(NB,:),'b+-','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
grid on;
xlabel('time(s)')
ylabel('Displacement(m)')
set(gca, 'FontName','Times New Roman','FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
set(get(gca,'ZLabel'),'Fontsize',20)
title('Dynamic displacement of note 3500');
legend('Displacement in x','Displacement in y')
subplot(3,1,2)
plot(Dt:Dt:Dt*Ns*Cn,x(NB-50,:),'r+-','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
plot(Dt:Dt:Dt*Ns*Cn,y(NB-50,:),'b+-','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
grid on;
xlabel('time(s)')
ylabel('Displacement(m)')
set(gca, 'FontName','Times New Roman','FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
set(get(gca,'ZLabel'),'Fontsize',20)
title('Dynamic displacement of note 3250');
legend('Displacement in x','Displacement in y')
subplot(3,1,3)
plot(Dt:Dt:Dt*Ns*Cn,x(NB-100,:),'r+-','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
plot(Dt:Dt:Dt*Ns*Cn,y(NB-100,:),'b+-','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
grid on;
xlabel('time(s)')
ylabel('Displacement(m)')
set(gca, 'FontName','Times New Roman','FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
set(get(gca,'ZLabel'),'Fontsize',20)
title('Dynamic displacement of note 3000');
legend('Displacement in x','Displacement in y')