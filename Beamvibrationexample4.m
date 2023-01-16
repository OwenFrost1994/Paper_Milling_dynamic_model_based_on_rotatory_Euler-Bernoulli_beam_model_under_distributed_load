clear;
clc;
%%梁模型在理想切削力作用下的振动情况%%
%%首先生成梁模型的基本参数，计算矩阵%%
%%定截面悬臂梁的振动算例%%
DL=0.02;%%梁单元的长度,mm
E=200000000000;%%杨氏模量N/m^2
p=7800;%%密度kg/(m^3)
L0=70;%%梁总悬长，mm，0号节点不用管，值已知
L0=L0/1000;%%梁总悬长，m
L1=40;%%第一段截面的长度，mm
L1=L1/1000;%%第一段截面的长度，m
L2=L0-L1;%%第二段截面的长度，mm
D1=10;%%第一段直径，mm
D1=D1/1000;%%第一段直径，m
D2=8;%%第二段直径，mm
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
% figure(1)
% plot(l,I,'r-');
% grid on;
% hold on;

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
% figure(2)
% plot(l,S,'ro');
% grid on;
% hold on;

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
% figure(3)
% subplot(3,1,1)
% plot(l,ID(:,1),'ro');
% hold on;
% subplot(3,1,2)
% plot(l,ID(:,2),'bo');
% hold on;
% subplot(3,1,3)
% plot(l,ID(:,3),'go');
% hold on;
% grid on;

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
    K(i,:)=E*DL*(ID(i,1)*A1(i,:)/(12*DL^2)+ID(i,2)*A2(i,:)/(2*DL^3)+ID(i,3)*A3(i,:)/(DL^4));
end
%%获得质量矩阵%%
for i=1:1:NB
    M(i,i)=p*DL*S(i);
end

%%其次计算切削力
%%不考虑刀具振动和偏心%%

%%刀具几何参数%%
D=1000*D1;%%刀具半径
N=4;%%刀具齿数
B=pi/6;%%刀具螺旋角
Cp=2*pi/N;%%齿间角
%%微元长度即为上面梁单元的微元长度

%%材料参数%%材料AL7075%%
Ktc=951.751;%%切向剪切力系数
Kte=14.0371;%%切向刃口力系数
Krc=608.561;%%径向剪切力系数
Kre=16.5002;%%径向刃口力系数
Kac=288.478;%%轴向剪切力系数
Kae=-1.25118;%%轴向刃口力系数

%%加工参数%%
%%铣削方式:顺铣%%
Cm=1;%%铣削方式，顺铣为1，逆铣为0
S=2000;%%主轴转速
f=600;%%进给速度
fs=10000;%%采样频率
ap=10;%%轴向切深（单位mm）
ae=3;%%径向切深（单位mm）
Cn=1;%%圈数circle number

%%基本参数计算%%
R=D/2;%%刀具半径
kb=(2*tan(B))/D;%%kβ计算
fe=f/(N*S);%%feed every tooth
w=2*pi*S/60;%%刀具角速度
T=2*pi/w;%%刀具周期
Ns=floor(60*fs/S);%%一个周期内的采样点个数
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
FY=zeros(NB,Ns*Cn);
FZ=zeros(NB,Ns*Cn);%%记录每一个梁单元上每一时刻的受力，其实际上只有头部的梁单元受力
Fx=0;
Fy=0;
Fz=0;
F=zeros(3,Ns*Cn);%%存储三个方向的切削力
x=zeros(NB,Ns*Cn);
y=zeros(NB,Ns*Cn);
K1=(M/(Dt^2)+K);
K1=inv(K1);
%%计算整个模拟域的大循环
%%计算切削力同步用有限差分法计算梁的振动位移和速度%%
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
            fx=(-cos(C))*(Ktc*fa+Kte)*Dl+(-sin(C))*(Krc*fa+Kre)*Dl;
            fy=( sin(C))*(Ktc*fa+Kte)*Dl+(-cos(C))*(Krc*fa+Kre)*Dl;
            fz=(Kac*fa+Kae)*Dl;
            else
                fx=0;fy=0;fz=0;
            end
            Fx=Fx+fx;Fy=Fy+fy;Fz=Fz+fz;
        end
        FX(NB-m+1,i)=Fx;FY(NB-m+1,i)=Fy;FZ(NB-m+1,i)=Fz;
        Fx=0;Fy=0;Fz=0;%%刀具受力累加变量归零
    end
    F(1,i)=sum(FX(:,i));F(2,i)=sum(FY(:,i));F(3,i)=sum(FZ(:,i));%%在矩阵中存储切削力
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

%%绘制力的图像，蓝色X向，红色Y向，绿色Z向
figure(5)
plot(Dt:Dt:Cn*Ns*Dt,F(1,:),'bo');
hold on;
plot(Dt:Dt:Cn*Ns*Dt,F(2,:),'ro');
hold on;
plot(Dt:Dt:Cn*Ns*Dt,F(3,:),'go');
grid on;

figure(6)
surf(x)
grid on;
hold on;
figure(7)
surf(y)
grid on;
hold on;