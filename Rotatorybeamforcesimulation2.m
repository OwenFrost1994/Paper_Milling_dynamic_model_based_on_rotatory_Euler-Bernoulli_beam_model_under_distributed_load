clear;
clc;
close all;
%%梁模型在理想切削力作用下的振动情况：切削力仿真+刀具弯曲情况%%
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


%%刀具几何参数%%
D=1000*D1;%%刀具半径
N=2;%%刀具齿数
B=pi*35/180;%%刀具螺旋角
Cp=2*pi/N;%%齿间角
%%微元长度即为上面梁单元的微元长度

%%加工参数%%
%%铣削方式:顺铣%%
Cm=1;%%铣削方式，顺铣为1，逆铣为0
SS=1000;%%主轴转速
f=80;%%进给速度
fs=20000;%%采样频率
ap=3;%%轴向切深（单位mm）
ae=2;%%径向切深（单位mm）
Cn=10;%%圈数circle number

Fm=xlsread('AE23.xlsx',4);
% Fm=downsample(Fm,20);
[m,n]=size(Fm);
A=zeros(m,1);
A=Fm(:,2);
Fm(:,2)=Fm(:,3);
Fm(:,3)=A;%%xy方向交换
clear A;
% Fm(:,3)=-Fm(:,3);
% Fm(:,2)=Fm(:,2)-1;
% Fm(:,3)=Fm(:,3)-2;
figure(2)
subplot(1,2,1)
plot(Fm(:,1),Fm(:,2),'b-','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
grid on;
subplot(1,2,2)
plot(Fm(:,1),Fm(:,3),'r-','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
grid on;
xlim([min(Fm(:,1)) max(Fm(:,1))])
xlabel('time(s)')
ylabel('Force(N)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('Actual forces in process');

%%材料参数%%材料AL%%
% Ktc=-67.802*ae^3+639.202*ae^2-1716.668*ae+3689.507;%%切向剪切力系数
% Kte=2.087*ae^3-16.463*ae^2+28.993*ae+11.804;%%切向刃口力系数
% Krc=61.124*ae^3-528.258*ae^2+1210.284*ae+337.690;%%径向剪切力系数
% Kre=-0.698*ae^3+10.441*ae^2-41.471*ae+39.937;%%径向刃口力系数
% Ktc=2684.35;%%ae=6
% Kte=18.07;
% Krc=2763.57;
% Kre=-20.91;

% Ktc=2461.67;%%ae=3
% Kte=6.96;
% Krc=864.56;
% Kre=-9.34;

Ktc=1161.15;%% ae=1 1632.47
Kte=7.16;
Krc=330.50;
Kre=-2.97;  

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

%%各种存储单元%%
Dl=1000*DL;
FX=zeros(NB,Ns*Cn);
FY=zeros(NB,Ns*Cn);%%记录每一个梁单元上每一时刻的受力，其实际上只有头部的梁单元受力
FXs=zeros(NB,Ns*Cn);
FYs=zeros(NB,Ns*Cn);
Fx=0;
Fy=0;
Fxs=0;
Fys=0;
F=zeros(2,Ns*Cn);%%存储两个方向的切削力
Fs=zeros(2,Ns*Cn);%%存储两个方向的静态切削力
x=zeros(NB,Ns*Cn);
y=zeros(NB,Ns*Cn);
Dx=zeros(ap/Dl,Ns*Cn);%%刀具在x方向的位移变化量，fv的第一个参数,这里动态切屑厚度只考虑切削段的部分，编号方式和刃编号方式一致
Dy=zeros(ap/Dl,Ns*Cn);%%刀具在y方向的位移变化量，fv的第二个参数
apx=zeros(1,Cn*N);
apy=zeros(1,Cn*N);
K1=(M/(Dt^2)+KS);
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
            fa=fe*sin(C)+Dx(m,i)*sin(C)+Dy(m,i)*cos(C);%未变形切屑厚度
            if C<Cex&&C>Cst&&fa>=0
            fx=(-cos(C))*(Ktc*fa+Kte)*Dl+(-sin(C))*(Krc*fa+Kre)*Dl;
            fy=( sin(C))*(Ktc*fa+Kte)*Dl+(-cos(C))*(Krc*fa+Kre)*Dl;
            else
                fx=0;fy=0;
            end
            Fx=Fx+fx;Fy=Fy+fy;
        end
        FX(NB-m+1,i)=Fx;FY(NB-m+1,i)=Fy;
        Fx=0;Fy=0;%%刀具受力累加变量归零
    end
    F(1,i)=sum(FX(:,i));F(2,i)=sum(FY(:,i));%%在矩阵中存储切削力
    %%x方向的振动
    if  i==1 %%判断是否是初始点（第一个点）
        F1=FX(:,i);
    else if i==2%%不是初始点
            F1=FX(:,i)-M*(0-2*x(:,i-1))/(Dt^2);
        else
            F1=FX(:,i)-M*(x(:,i-2)-2*x(:,i-1))/(Dt^2);
        end
    end
    x(:,i)=K1*F1;
    %%y方向的振动
    if  i==1 %%判断是否是初始点（第一个点）
        F2=FY(:,i);
    else if i==2%%不是初始点
            F2=FY(:,i)-M*(0-2*y(:,i-1))/(Dt^2);
        else
            F2=FY(:,i)-M*(y(:,i-2)-2*y(:,i-1))/(Dt^2);
        end
    end
    y(:,i)=K1*F2;
end
Fd=F;

Ca=Cs;%%初始角度
F=zeros(2,Ns*Cn);
Dx=zeros(ap/Dl,Ns*Cn);
Dy=zeros(ap/Dl,Ns*Cn);
apx=zeros(1,Cn*N);
apy=zeros(1,Cn*N);
x=zeros(NB,Ns*Cn);
y=zeros(NB,Ns*Cn);

%%计算广义动刚度矩阵和广义阻尼矩阵，建立广义刚度和广义质量矩阵
[i,j]=size(KS);
KW=-w^2*p*diag([S,S]);
CQ=[zeros(i,j),diag(-2*w*p*S*DL/1000);diag(2*w*p*S*DL/1000),zeros(i,j)];
KQ=[KS,zeros(i,j);zeros(i,j),KS]+KW;
MQ=[M,zeros(i,j);zeros(i,j),M];
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
        else
            fx=0;fy=0;
        end
        Fx=Fx+fx;Fy=Fy+fy;%%累加每一个刀片上每个齿所受的切削力
        end
        FX(NB-m+1,i)=Fx;FY(NB-m+1,i)=Fy;
        Fx=0;Fy=0;%%刀具受力累加变量归零
    end
    F(1,i)=sum(FX(:,i));F(2,i)=sum(FY(:,i));%%在矩阵中存储切削力
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
FR=F;

for i=1:1:Ns*Cn;
    Ca=Ca+DC;%%微元角度叠加计算刀具的转动角
    if Ca>=2*pi%%考虑刀具多个旋转周期，累加的刀具角度超过一周就减去一个2π
        Ca=Ca-2*pi;
    else
    end
    
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
            fa=fe*sin(C);%未变形切屑厚度
            if C<Cex&&C>Cst&&fa>=0
            fxs=(-cos(C))*(Ktc*fa+Kte)*Dl+(-sin(C))*(Krc*fa+Kre)*Dl;
            fys=( sin(C))*(Ktc*fa+Kte)*Dl+(-cos(C))*(Krc*fa+Kre)*Dl;
            else
                fxs=0;fys=0;
            end
            Fxs=Fxs+fxs;Fys=Fys+fys;
        end
        FXs(NB-m+1,i)=Fxs;FYs(NB-m+1,i)=Fys;
        Fxs=0;Fys=0;%%刀具受力累加变量归零
    end
    Fs(1,i)=sum(FXs(:,i));Fs(2,i)=sum(FYs(:,i));%%在矩阵中存储切削力
end
%%绘制力的图像，蓝色X向，红色Y向，绿色Z向
% figure(2)
% plot(Dt:Dt:Cn*Ns*Dt,F(1,:),'bo');
% hold on;
% plot(Dt:Dt:Cn*Ns*Dt,F(2,:),'ro');
% hold on;

%%力的对比加进来%%
n1=5*Ns;
n2=9*Ns;
n3=4*Ns-200;
n4=8*Ns-200;
figure(3)
% title('Measured and predicted cutting forces','Fontsize',20)
subplot(2,1,1)
plot(Fm(n3:n4,1),FR(1,n1:n2),'m-','Markersize',7,'Markerface','white','linewidth',3.0);%%这里要注意F是3×N,读取的数据是N×4
hold on;
plot(Fm(n3:n4,1),Fd(1,n1:n2),'r-','Markersize',7,'Markerface','white','linewidth',1.0);%%这里要注意F是3×N,读取的数据是N×4
hold on;
plot(Fm(n3:n4,1),Fs(1,n1:n2),'b-.','Markersize',7,'Markerface','white','linewidth',1.0);
hold on;
plot(Fm(n3:n4,1),Fm(n3:n4,2),'k--','Markersize',7,'Markerface','white','linewidth',1.0);
hold on;
grid on;
xlim([min(Fm(n3:n4,1)) max(Fm(n3:n4,1))])
xlabel('time(s)')
ylabel('Fx(N)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
legend('Rorarory model','Static model','Traditional method','Measured signal')
subplot(2,1,2)
plot(Fm(n3:n4,1),FR(2,n1:n2),'m-','Markersize',7,'Markerface','white','linewidth',3.0);%%这里要注意F是3×N,读取的数据是N×4
hold on;
plot(Fm(n3:n4,1),Fd(2,n1:n2),'r-','Markersize',7,'Markerface','white','linewidth',1.0);
hold on;
plot(Fm(n3:n4,1),Fs(2,n1:n2),'b-.','Markersize',7,'Markerface','white','linewidth',1.0);
hold on;
plot(Fm(n3:n4,1),Fm(n3:n4,3),'k--','Markersize',7,'Markerface','white','linewidth',1.0);
hold on;
grid on;
xlim([min(Fm(n3:n4,1)) max(Fm(n3:n4,1))])
xlabel('time(s)')
ylabel('Fy(N)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
legend('Rorarory model','Static model','Traditional method','Measured signal')

Fb=Fm(n3:n4,2:4);
Fb(1:Ns/2,:)=(Fb(1:Ns/2,:)+Fb(Ns/2+1:Ns,:))/2;
Fb(Ns/2+1:Ns,:)=Fb(1:Ns/2,:);
Fb(Ns+1:Ns+Ns/2,:)=(Fb(Ns+1:Ns+Ns/2,:)+Fb(Ns+Ns/2+1:2*Ns,:))/2;
Fb(Ns+Ns/2+1:2*Ns,:)=Fb(Ns+1:Ns+Ns/2,:);
Fb(2*Ns+1:2*Ns+Ns/2,:)=(Fb(2*Ns+1:2*Ns+Ns/2,:)+Fb(2*Ns+Ns/2+1:3*Ns,:))/2;
Fb(2*Ns+Ns/2+1:3*Ns,:)=Fb(2*Ns+1:2*Ns+Ns/2,:);
Fb(3*Ns+1:3*Ns+Ns/2,:)=(Fb(3*Ns+1:3*Ns+Ns/2,:)+Fb(3*Ns+Ns/2+1:4*Ns,:))/2;
Fb(3*Ns+Ns/2+1:4*Ns,:)=Fb(3*Ns+1:3*Ns+Ns/2,:);

figure(4)
% title('Measured and predicted cutting forces','Fontsize',20)
subplot(2,1,1)
plot(Fm(n3:n4,1),FR(1,n1:n2),'m-','Markersize',7,'Markerface','white','linewidth',3.0);%%这里要注意F是3×N,读取的数据是N×4
hold on;
plot(Fm(n3:n4,1),Fd(1,n1:n2),'r-','Markersize',7,'Markerface','white','linewidth',1.0);
hold on;
plot(Fm(n3:n4,1),Fs(1,n1:n2),'b-.','Markersize',7,'Markerface','white','linewidth',1.0);
hold on;
plot(Fm(n3:n4,1),Fb(:,1),'k--','Markersize',7,'Markerface','white','linewidth',1.0);
hold on;
grid on;
xlim([min(Fm(n3:n4,1)) max(Fm(n3:n4,1))])
xlabel('time(s)')
ylabel('Fx(N)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
legend('Rorarory model','Static model','Traditional method','Balanced signal')
subplot(2,1,2)
plot(Fm(n3:n4,1),FR(2,n1:n2),'m-','Markersize',7,'Markerface','white','linewidth',3.0);%%这里要注意F是3×N,读取的数据是N×4
hold on;
plot(Fm(n3:n4,1),Fd(2,n1:n2),'r-','Markersize',7,'Markerface','white','linewidth',1.0);
hold on;
plot(Fm(n3:n4,1),Fs(2,n1:n2),'b-.','Markersize',7,'Markerface','white','linewidth',1.0);
hold on;
plot(Fm(n3:n4,1),Fb(:,2),'k--','Markersize',7,'Markerface','white','linewidth',1.0);
hold on;
grid on;
xlim([min(Fm(n3:n4,1)) max(Fm(n3:n4,1))])
xlabel('time(s)')
ylabel('Fy(N)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
legend('Rorarory model','Static model','Traditional method','Balanced signal')

