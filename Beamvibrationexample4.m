clear;
clc;
%%��ģ�������������������µ������%%
%%����������ģ�͵Ļ����������������%%
%%��������������������%%
DL=0.02;%%����Ԫ�ĳ���,mm
E=200000000000;%%����ģ��N/m^2
p=7800;%%�ܶ�kg/(m^3)
L0=70;%%����������mm��0�Žڵ㲻�ùܣ�ֵ��֪
L0=L0/1000;%%����������m
L1=40;%%��һ�ν���ĳ��ȣ�mm
L1=L1/1000;%%��һ�ν���ĳ��ȣ�m
L2=L0-L1;%%�ڶ��ν���ĳ��ȣ�mm
D1=10;%%��һ��ֱ����mm
D1=D1/1000;%%��һ��ֱ����m
D2=8;%%�ڶ���ֱ����mm
D2=D2/1000;%%�ڶ���ֱ����m
I1=pi*(D1)^4/64;%%��һ�ι��Ծ�,m^4
I2=pi*(D2)^4/64;%%�ڶ��ι��Ծ�,m^4
DL=DL/1000;%%����Ԫ�ĳ���,m
NB=ceil(L0/DL);%%����Ԫ�ĸ���

%%���Ծ���ɢ%%
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

%%�������ɢ%%
%%�����i�Ǵ�1��Nһ��N���ڵ㣬û����ԭ�нڵ�Ļ����ϼ�����ڵ�
for i=1:1:NB
    l(i)=i*DL-0.5*DL;%%�ڵ��λ�ô���ÿ��Ƭ����е�
    if l(i)<=L1%%ǰ��ν����
        S(i)=pi*D1^2/4;
    else
        S(i)=pi*D2^2/4;%%���ν����
    end
end
% figure(2)
% plot(l,S,'ro');
% grid on;
% hold on;

%%����I'',I',I������ID�У�ÿһ���ڵ���������ֵ����������
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

%%����Ϊ׼������%%
%%�������ÿ��ʱ������ϵ�����󣬳������󣬺;������%%
%%����ǰ������ڵ��Ѿ��ڹ�ʽ���п��ǹ���%%
%%����ϵ������%%
A1=zeros(NB,NB);
A2=zeros(NB,NB);
A3=zeros(NB,NB);
%%���׵���ϵ������%%
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
%%���׵���ϵ������%%
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
%%�Ľ׵���ϵ������%%
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
%%��øնȾ���%%
for i=1:1:NB
    K(i,:)=E*DL*(ID(i,1)*A1(i,:)/(12*DL^2)+ID(i,2)*A2(i,:)/(2*DL^3)+ID(i,3)*A3(i,:)/(DL^4));
end
%%�����������%%
for i=1:1:NB
    M(i,i)=p*DL*S(i);
end

%%��μ���������
%%�����ǵ����񶯺�ƫ��%%

%%���߼��β���%%
D=1000*D1;%%���߰뾶
N=4;%%���߳���
B=pi/6;%%����������
Cp=2*pi/N;%%�ݼ��
%%΢Ԫ���ȼ�Ϊ��������Ԫ��΢Ԫ����

%%���ϲ���%%����AL7075%%
Ktc=951.751;%%���������ϵ��
Kte=14.0371;%%�����п���ϵ��
Krc=608.561;%%���������ϵ��
Kre=16.5002;%%�����п���ϵ��
Kac=288.478;%%���������ϵ��
Kae=-1.25118;%%�����п���ϵ��

%%�ӹ�����%%
%%ϳ����ʽ:˳ϳ%%
Cm=1;%%ϳ����ʽ��˳ϳΪ1����ϳΪ0
S=2000;%%����ת��
f=600;%%�����ٶ�
fs=10000;%%����Ƶ��
ap=10;%%���������λmm��
ae=3;%%���������λmm��
Cn=1;%%Ȧ��circle number

%%������������%%
R=D/2;%%���߰뾶
kb=(2*tan(B))/D;%%k�¼���
fe=f/(N*S);%%feed every tooth
w=2*pi*S/60;%%���߽��ٶ�
T=2*pi/w;%%��������
Ns=floor(60*fs/S);%%һ�������ڵĲ��������
if Cm==1%%˳ϳ
    Cst=pi-acos((R-ae)/R);%%�����
    Cex=pi;%%�г���
else%%��ϳ
    Cst=0;%%�����
    Cex=acos((R-ae)/R);%%�г���
end
Cs=0;%%��ʼ�Ƕ�
Dt=T/Ns;%%ʱ�䲽��
DC=Dt*w;%%�Ƕ�����
Ca=Cs;%%��ʼ�Ƕ�

%%���ִ洢��Ԫ%%
Dl=1000*DL;
FX=zeros(NB,Ns*Cn);
FY=zeros(NB,Ns*Cn);
FZ=zeros(NB,Ns*Cn);%%��¼ÿһ������Ԫ��ÿһʱ�̵���������ʵ����ֻ��ͷ��������Ԫ����
Fx=0;
Fy=0;
Fz=0;
F=zeros(3,Ns*Cn);%%�洢���������������
x=zeros(NB,Ns*Cn);
y=zeros(NB,Ns*Cn);
K1=(M/(Dt^2)+K);
K1=inv(K1);
%%��������ģ����Ĵ�ѭ��
%%����������ͬ�������޲�ַ�����������λ�ƺ��ٶ�%%
for i=1:1:Ns*Cn;
    Ca=Ca+DC;%%΢Ԫ�Ƕȵ��Ӽ��㵶�ߵ�ת����
    if Ca>=2*pi%%���ǵ��߶����ת���ڣ��ۼӵĵ��߽Ƕȳ���һ�ܾͼ�ȥһ��2��
        Ca=Ca-2*pi;
    else
    end
    
    %%΢Ԫ���Ӽ���һ�������ϵ�������
    for m=1:1:ap/Dl
        Cd=Ca-m*kb*Dl+0.5*kb*Dl;
        %���Ӽ��������ݵ���ѭ��
        for j=1:1:N
        C=Cd-(j-1)*Cp;%%���Ƕ�ݴ��ڵĳݼ���ͺ�
        if C<0
            C=C+2*pi;%%������ݽǶ�С������ת�����ĽǶ�
        else
        end

            fa=fe*sin(C);%δ������м���
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
        Fx=0;Fy=0;Fz=0;%%���������ۼӱ�������
    end
    F(1,i)=sum(FX(:,i));F(2,i)=sum(FY(:,i));F(3,i)=sum(FZ(:,i));%%�ھ����д洢������
    %%x�������
        if  i==1 %%�ж��Ƿ��ǳ�ʼ�㣨��һ���㣩
        F1=FX(:,i);
    else if i==2%%���ǳ�ʼ��
            F1=FX(:,i)-M*(0-2*x(:,i-1))/(Dt^2);
        else
            F1=FX(:,i)-M*(x(:,i-2)-2*x(:,i-1))/(Dt^2);
        end
    end
    x(:,i)=K1*F1;
    %%y�������
    if  i==1 %%�ж��Ƿ��ǳ�ʼ�㣨��һ���㣩
        F2=FY(:,i);
    else if i==2%%���ǳ�ʼ��
            F2=FY(:,i)-M*(0-2*y(:,i-1))/(Dt^2);
        else
            F2=FY(:,i)-M*(y(:,i-2)-2*y(:,i-1))/(Dt^2);
        end
    end
    y(:,i)=K1*F2;
end

%%��������ͼ����ɫX�򣬺�ɫY����ɫZ��
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