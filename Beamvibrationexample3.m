clear;
clc;
%%变截面悬臂梁的振动算例%%
%%考虑到头尾节点的边界条件这里的梁实际上包括头部和尾部的虚拟节点%%
DL=0.5;%%梁单元的长度,mm
E=200000000000;%%杨氏模量N/m^2
p=7800;%%密度kg/(m^3)
L0=70+1*DL;%%梁总悬长，mm，0号节点不用管，值已知
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
N=ceil(L0/DL);%%梁单元的个数

%%惯性矩离散%%
%%这里的i是从1到N+1一共，N+1个节点，也就是在原有节点的基础上加了一个虚拟节点
for i=1:1:N
    l(i)=i*DL-0.5*DL;
    if l(i)<=L1%%前半段惯性矩
        I(i)=I1;
    else
        I(i)=I2;%%后半段惯性矩
    end
end
% figure(1)
% plot(l,I,'ro');
% grid on;
% hold on;

%%截面积离散%%
%%这里的i是从1到N+1一共，N+1个节点，也就是在原有节点的基础上加了一个虚拟节点
for i=1:1:N
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

%%梁单元受力情况
L=10;%%受力部分的长度，mm
L=L/1000;%%受力部分的长度，m
Tw=2;%%计算总时域，s
DT=0.01;%%时间步长，s
M=Tw/DT;%%时间节点个数
%%受到的力为时变力
for j=1:1:M
    t(j)=DT*j;
    for i=1:1:N
        if l(i)<=L0-L-DL
            F(i,j)=0;%%一部分长度上受力为0，包括最后一个节点
        else if l(i)<=L0-DL
                F(i,j)=30*(L0-l(i))/L*DL*sin(2*pi*t(j));
            else F(i,j)=0;
            end
        end
    end
end
figure(3)
surf(F)
grid on;
hold on;
%%计算I'',I',I，存入ID中，每一个节点的这个三个值都是行向量
for i=1:1:N
    if i==1
        I2=(I(i)-2*I(i+1)+I(i+2))/(DL^2);
        I1=(-3*I(i)+4*I(i+1)-I(i+2))/(2*DL);
    else if i==N
                    I2=(I(i-2)-2*I(i-1)+I(i))/(DL^2);
                    I1=(I(i-2)-4*I(i-1)+3*I(i))/(2*DL);
        else
            I2=(I(i-1)-2*I(i)+I(i+1))/(DL^2);
            I1=(I(i+1)-I(i-1))/(2*DL);
        end
    end
    ID(i,:)=[I2,2*I1,I(i)];
end
% figure(4)
% subplot(3,1,1)
% plot(l,EID(:,1),'ro');
% hold on;
% subplot(3,1,2)
% plot(l,EID(:,2),'bo');
% hold on;
% subplot(3,1,3)
% plot(l,EID(:,3),'go');
% hold on;
% grid on;

%%以上为准备工作%%
%%下面根据每个时刻生成系数矩阵，常数矩阵，和矩阵求解%%
%%这里前后虚拟节点已经在公式当中考虑过了%%
%%产生系数矩阵%%
A1=zeros(N,N);
A2=zeros(N,N);
A3=zeros(N,N);
%%二阶导的系数矩阵%%
for i=1:1:N
    if i==1
        A1(i,1:5)=[-31,16,-1,0,0];
    else if i==2
            A1(i,1:5)=[16,-30,16,-1,0];
        else if i==N-1
                A1(i,N-4:N)=[0,-1,111/7,-201/7,97/7];
            else if i==N
                    A1(i,N-4:N)=[0,0,0,0,0];
                else
                    A1(i,i-2:i+2)=[-1,16,-30,16,-1];
                end
            end
        end
    end
end
%%三阶导的系数矩阵%%
for i=1:1:N
    if i==1
        A2(i,1:5)=[-1,-2,1,0,0];
    else if i==2
            A2(i,1:5)=[2,0,-2,1,0];
        else if i==N-1
                A2(i,N-4:N)=[0,-1,15/7,-9/7,1/7];
            else if i==N
                    A2(i,N-4:N)=[0,0,0,0,0];
                else
                    A2(i,i-2:i+2)=[-1,2,0,2,-1];
                end
            end
        end
    end
end
%%四阶导的系数矩阵%%
for i=1:1:N
    if i==1
        A3(i,1:5)=[7,-4,1,0,0];
    else if i==2
            A3(i,1:5)=[-4,6,-4,1,0];
        else if i==N-1
                A3(i,N-4:N)=[0,1,-27/7,33/7,-13/7];
            else if i==N
                    A3(i,N-4:N)=[0,0,12/7,-24/7,12/7];
                else
                    A3(i,i-2:i+2)=[1,-4,6,-4,1];
                end
            end
        end
    end
end
%%获得刚度矩阵%%
for i=1:1:N
    K(i,:)=E*DL*(ID(i,1)*A1(i,:)/(12*DL^2)+ID(i,2)*A2(i,:)/(2*DL^3)+ID(i,3)*A3(i,:)/(DL^4));
end
%%获得质量矩阵%%
for i=1:1:N
    M(i,i)=p*DL*S(i);
end
Mt=Tw/DT;%%时间节点个数
y=zeros(N,Mt);
y1=zeros(N,Mt);
y0=zeros(N,1);
y10=zeros(N,1);
%%Runge-Kutta法计算梁的振动位移和速度%%
M1=inv(M);
M2=M1*K;
for j=1:1:Mt
    if  j==1 %%判断是否是初始点（第一个点）
        K1=y10;L1=-M\K*y0+M\zeros(N,1);
        K2=y10+DT*L1/2;L2=-M\K*(y0+DT*K1/2)+M\((zeros(N,1)+F(:,j))/2);
        K3=y10+DT*L2/2;L3=-M\K*(y0+DT*K2/2)+M\((zeros(N,1)+F(:,j))/2);
        K4=y10+DT*L3;L4=-M\K*(y0+DT*K3)+M\F(:,j);
        y(:,j)=y0+DT*(K1+2*K2+2*K3+K4)/6;
        y1(:,j)=y10+DT*(L1+2*L2+2*L3+L4)/6;
    else %%不是初始点
        K1=y1(:,j-1);L1=-M\K*y(:,j-1)+M\F(:,j-1);
        K2=y1(:,j-1)+DT*L1/2;L2=-M\K*(y(:,j-1)+DT*K1/2)+M\((F(:,j-1)+F(:,j))/2);
        K3=y1(:,j-1)+DT*L2/2;L3=-M\K*(y(:,j-1)+DT*K2/2)+M\((F(:,j-1)+F(:,j))/2);
        K4=y1(:,j-1)+DT*L3;L4=-M\K*(y(:,j-1)+DT*K3)+M\F(:,j);
        y(:,j)=y(:,j-1)+DT*(K1+2*K2+2*K3+K4)/6;
        y1(:,j)=y(:,j-1)+DT*(L1+2*L2+2*L3+L4)/6;
    end
end
% figure(5)
% surf(y)
% grid on;
% hold on;