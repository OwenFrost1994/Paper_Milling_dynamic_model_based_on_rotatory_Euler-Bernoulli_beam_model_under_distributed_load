clear;
clc;
%%变截面悬臂梁的振动算例%%
E=140000000000;%%杨氏模量
p=7800;%%密度kg/(m^3)
L0=70;%%梁总悬长，mm
L0=L0/1000;%%梁总悬长，m
L1=40;%%第一段截面的长度，mm
L1=L1/1000;%%第一段截面的长度，m
L2=L0-L1;%%第二段截面的长度，m
D1=10;%%第一段直径，mm
D1=D1/1000;%%第一段直径，m
D2=8;%%第二段直径，mm
D2=D2/1000;%%第二段直径，m
I1=pi*(D1)^4/64;%%第一段惯性矩
I2=pi*(D2)^4/64;%%第二段惯性矩
DL=1;%%梁单元的长度,mm
DL=DL/1000;%%梁单元的长度,m
N=L0/DL;%%梁单元的个数

%%惯性矩离散%%
for i=1:1:N
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
for i=1:1:N
    l(i)=i*DL-0.5*DL;
    if l(i)<=L1
        S(i)=pi*D1^2/4;
    else
        S(i)=pi*D2^2/4;
    end
end
figure(2)
plot(l,S,'r-');
grid on;
hold on;

%%梁单元受力情况
L=20;%%受力部分的长度，mm
L=L/1000;%%受力部分的长度，m
Tw=2;%%计算总时域，s
DT=0.01;%%时间步长，s
M=Tw/DT;%%时间节点个数
%%受到的力为时变力
for j=1:1:M
    t(j)=DT*j;
    for i=1:1:N
        if l(i)<=L0-L
            F(i,j)=0;%%一部分长度上受力为0
        else
            F(i,j)=30000*(L0-l(i))/L*DL*sin(2*pi*t(j));
        end
    end
end
figure(3)
surf(F)
grid on;
hold on;
%%计算EI'',EI',EI，存入EID中，每一个节点的这个三个值都是行向量
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
    EID(i,:)=E*[I2,2*I1,I(i)];
end
figure(4)
subplot(3,1,1)
plot(l,EID(:,1),'ro');
hold on;
subplot(3,1,2)
plot(l,EID(:,2),'bo');
hold on;
subplot(3,1,3)
plot(l,EID(:,3),'go');
hold on;
grid on;

%%以上为准备工作%%
%%下面根据每个时刻生成系数矩阵，常数矩阵，和矩阵求解%%

%%先假定某一时刻t(j)
j=20;
A=zeros(N,N);
y=zeros(N,M);
%%这里没有加边值条件，首先肯定不收敛%%
for i=1:1:N
    if i==1%%头节点，i=1
        A(i,1:5)=EID(i,:)*[35/12,-26/3,19/2,-14/3,11/12;-5/2,9,-12,7,-3/2;1,-4,6,-4,1]+[p*S(i)/(3*DT^2),0,0,0,0];
    else if i==2%%头节点，i=2
            A(i,1:5)=EID(i,:)*[11/12,-5/3,1/2,1/3,-1/12;-3/2,5,-6,3,-1/2;1,-4,6,-4,1]+[0,p*S(i)/(3*DT^2),0,0,0];
        else if i==N%%尾节点，i=N
                A(i,N-4:N)=EID(i,:)*[11/12,-14/3,19/2,-26/3,35/12;2/3,-7,12,-9,10/3;1,-4,6,-4,1]+[0,0,0,0,p*S(i)/(3*DT^2)];
            else if i==N-1%%尾节点，i=N-1
                    A(i,N-4:N)=EID(i,:)*[-1/12,1/3,1/2,-5/3,11/12;1/2,-3,6,-5,3/2;1,-4,6,-4,1]+[0,0,0,p*S(i)/(3*DT^2),0];
                else
                    A(i,i-2:i+2)=EID(i,:)*[-1/12,4/3,-5/2,4/3,-1/12;-1/2,1,0,-1,1/2;1,-4,6,-4,1]+[0,0,p*S(i)/(3*DT^2),0,0];%%普通中间节点i
                end
            end
        end
    end
end
for j=1:1:M
if j==1
    B=F(:,j)-0;
else if j==2
        B=F(:,j)-p*S(:).*(0-2*y(:,j-1))/(3*DT^2);
    else
        B=F(:,j)-p*S(:).*(y(:,j-2)-2*y(:,j-1))/(3*DT^2);
    end
end
y(:,j)=A\B;
end
figure(5)
surf(y)
grid on;
hold on;