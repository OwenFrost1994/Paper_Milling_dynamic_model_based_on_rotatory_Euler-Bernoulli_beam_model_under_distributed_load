clear;
clc;
close all;
%%��ģ�������������������µ������������������+�����������%%
%%����������ģ�͵Ļ����������������%%
%%���ߵ�ЧΪ���������ξ�̬��ģ��%%
DL=0.02;%%����Ԫ�ĳ���,mm
E=577500000000;%%����ģ��N/m^2
p=14500;%%�ܶ�kg/(m^3)
L0=60;%%����������mm��0�Žڵ㲻�ùܣ�ֵ��֪
L0=L0/1000;%%����������m
L1=30;%%��һ�ν���ĳ��ȣ�mm
L1=L1/1000;%%��һ�ν���ĳ��ȣ�m
L2=L0-L1;%%�ڶ��ν���ĳ��ȣ�mm
D1=6;%%��һ��ֱ����mm
D1=D1/1000;%%��һ��ֱ����m
D2=4.8;%%�ڶ���ֱ����mm
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
    KS(i,:)=E*DL*(ID(i,1)*A1(i,:)/(12*DL^2)+ID(i,2)*A2(i,:)/(2*DL^3)+ID(i,3)*A3(i,:)/(DL^4));
end
%%�����������%%
for i=1:1:NB
    M(i,i)=p*DL*S(i);
end


%%���߼��β���%%
D=1000*D1;%%���߰뾶
N=2;%%���߳���
B=pi*35/180;%%����������
Cp=2*pi/N;%%�ݼ��
%%΢Ԫ���ȼ�Ϊ��������Ԫ��΢Ԫ����

%%�ӹ�����%%
%%ϳ����ʽ:˳ϳ%%
Cm=1;%%ϳ����ʽ��˳ϳΪ1����ϳΪ0
SS=1000;%%����ת��
f=80;%%�����ٶ�
fs=20000;%%����Ƶ��
ap=3;%%���������λmm��
ae=2;%%���������λmm��
Cn=10;%%Ȧ��circle number

Fm=xlsread('AE23.xlsx',4);
% Fm=downsample(Fm,20);
[m,n]=size(Fm);
A=zeros(m,1);
A=Fm(:,2);
Fm(:,2)=Fm(:,3);
Fm(:,3)=A;%%xy���򽻻�
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

%%���ϲ���%%����AL%%
% Ktc=-67.802*ae^3+639.202*ae^2-1716.668*ae+3689.507;%%���������ϵ��
% Kte=2.087*ae^3-16.463*ae^2+28.993*ae+11.804;%%�����п���ϵ��
% Krc=61.124*ae^3-528.258*ae^2+1210.284*ae+337.690;%%���������ϵ��
% Kre=-0.698*ae^3+10.441*ae^2-41.471*ae+39.937;%%�����п���ϵ��
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

%%������������%%
R=D/2;%%���߰뾶
kb=(2*tan(B))/D;%%k�¼���
fe=f/(N*SS);%%feed every tooth
w=2*pi*SS/60;%%���߽��ٶ�
T=2*pi/w;%%��������
Ns=floor(60*fs/SS);%%һ�������ڵĲ��������
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
FY=zeros(NB,Ns*Cn);%%��¼ÿһ������Ԫ��ÿһʱ�̵���������ʵ����ֻ��ͷ��������Ԫ����
FXs=zeros(NB,Ns*Cn);
FYs=zeros(NB,Ns*Cn);
Fx=0;
Fy=0;
Fxs=0;
Fys=0;
F=zeros(2,Ns*Cn);%%�洢���������������
Fs=zeros(2,Ns*Cn);%%�洢��������ľ�̬������
x=zeros(NB,Ns*Cn);
y=zeros(NB,Ns*Cn);
Dx=zeros(ap/Dl,Ns*Cn);%%������x�����λ�Ʊ仯����fv�ĵ�һ������,���ﶯ̬��м���ֻ���������εĲ��֣���ŷ�ʽ���б�ŷ�ʽһ��
Dy=zeros(ap/Dl,Ns*Cn);%%������y�����λ�Ʊ仯����fv�ĵڶ�������
apx=zeros(1,Cn*N);
apy=zeros(1,Cn*N);
K1=(M/(Dt^2)+KS);
K1=inv(K1);
%%��������ģ����Ĵ�ѭ��
%%����������ͬ�������޲�ַ�����������λ�ƺ��ٶ�%%
for i=1:1:Ns*Cn;
    Ca=Ca+DC;%%΢Ԫ�Ƕȵ��Ӽ��㵶�ߵ�ת����
    if Ca>=2*pi%%���ǵ��߶����ת���ڣ��ۼӵĵ��߽Ƕȳ���һ�ܾͼ�ȥһ��2��
        Ca=Ca-2*pi;
    else
    end
     %%���Ƕ�̬��м���%%
    for m=1:1:ap/Dl
        if i<=(Ns/N)%%����ǵ�һ���������е�ʱ����м�ϱ���û��ǰһ�������г��Ĳ��Ʊ��棬��ʱ�൱�ڵ������񶯵����
            Dx(m,i)=1000*(x(NB-m+1,i)-0);%%ע�⵽����ѧ�����������е�λΪ�����Ƶ�λ������õ�λ�Ƶ�λΪm����ÿ�ݽ������ĵ�λ��mm
            Dy(m,i)=1000*(y(NB-m+1,i)-0);%%ע��ʱ����ÿһ�У���ͬ���б������
        else
            q=floor(i*N/Ns);
            if (i-q*Ns/N)==0
                q=q-1;%%���������ĳ�����ݽ���׶Σ�ʵ������Ҫ��ȥһ������
            else
            end
            for n=1:1:q%%ѭ��ȡ��ǰ��Ľڵ�
                apx(n)=x(NB-m+1,i-n*Ns/N);%%ע���ص��᷽���ÿ�����ݶ�Ҫ����
                apy(n)=y(NB-m+1,i-n*Ns/N);%%ע���ص��᷽���ÿ�����ݶ�Ҫ����
            end
            Dx(m,i)=1000*(x(i)-min(apx(1:q)));
            Dy(m,i)=1000*(y(i)-min(apy(1:q)));
        end
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
            fa=fe*sin(C)+Dx(m,i)*sin(C)+Dy(m,i)*cos(C);%δ������м���
            if C<Cex&&C>Cst&&fa>=0
            fx=(-cos(C))*(Ktc*fa+Kte)*Dl+(-sin(C))*(Krc*fa+Kre)*Dl;
            fy=( sin(C))*(Ktc*fa+Kte)*Dl+(-cos(C))*(Krc*fa+Kre)*Dl;
            else
                fx=0;fy=0;
            end
            Fx=Fx+fx;Fy=Fy+fy;
        end
        FX(NB-m+1,i)=Fx;FY(NB-m+1,i)=Fy;
        Fx=0;Fy=0;%%���������ۼӱ�������
    end
    F(1,i)=sum(FX(:,i));F(2,i)=sum(FY(:,i));%%�ھ����д洢������
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
Fd=F;

Ca=Cs;%%��ʼ�Ƕ�
F=zeros(2,Ns*Cn);
Dx=zeros(ap/Dl,Ns*Cn);
Dy=zeros(ap/Dl,Ns*Cn);
apx=zeros(1,Cn*N);
apy=zeros(1,Cn*N);
x=zeros(NB,Ns*Cn);
y=zeros(NB,Ns*Cn);

%%������嶯�նȾ���͹���������󣬽�������նȺ͹�����������
[i,j]=size(KS);
KW=-w^2*p*diag([S,S]);
CQ=[zeros(i,j),diag(-2*w*p*S*DL/1000);diag(2*w*p*S*DL/1000),zeros(i,j)];
KQ=[KS,zeros(i,j);zeros(i,j),KS]+KW;
MQ=[M,zeros(i,j);zeros(i,j),M];
%%΢�ַ�����ɢ�����
K1=MQ/(Dt^2)+3*CQ/(2*Dt)+KQ;
K2=2*MQ/(Dt^2)+2*CQ/Dt;
K3=-MQ/(Dt^2)-CQ/2*Dt;
K1=inv(K1);
%%��������ģ����Ĵ�ѭ��
%%����������ͬ�������޲�ַ�����������λ�ƺ��ٶ�%%
for i=1:1:Ns*Cn;
    Ca=Ca+DC;%%΢Ԫ�Ƕȵ��Ӽ��㵶�ߵ�ת����
    if Ca>=2*pi%%���ǵ��߶����ת���ڣ��ۼӵĵ��߽Ƕȳ���һ�ܾͼ�ȥһ��2��
        Ca=Ca-2*pi;
    else
    end
    %%���Ƕ�̬��м���%%
    for m=1:1:ap/Dl
        if i<=(Ns/N)%%����ǵ�һ���������е�ʱ����м�ϱ���û��ǰһ�������г��Ĳ��Ʊ��棬��ʱ�൱�ڵ������񶯵����
            Dx(m,i)=1000*(x(NB-m+1,i)-0);%%ע�⵽����ѧ�����������е�λΪ�����Ƶ�λ������õ�λ�Ƶ�λΪm����ÿ�ݽ������ĵ�λ��mm
            Dy(m,i)=1000*(y(NB-m+1,i)-0);%%ע��ʱ����ÿһ�У���ͬ���б������
        else
            q=floor(i*N/Ns);
            if (i-q*Ns/N)==0
                q=q-1;%%���������ĳ�����ݽ���׶Σ�ʵ������Ҫ��ȥһ������
            else
            end
            for n=1:1:q%%ѭ��ȡ��ǰ��Ľڵ�
                apx(n)=x(NB-m+1,i-n*Ns/N);%%ע���ص��᷽���ÿ�����ݶ�Ҫ����
                apy(n)=y(NB-m+1,i-n*Ns/N);%%ע���ص��᷽���ÿ�����ݶ�Ҫ����
            end
            Dx(m,i)=1000*(x(i)-min(apx(1:q)));
            Dy(m,i)=1000*(y(i)-min(apy(1:q)));
        end
    end

    %%���ս���̬��м��ȴ��������%%
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
        fa=fe*sin(C)++Dx(m,i)*sin(C)+Dy(m,i)*cos(C);%δ������м���
        if C<Cex&&C>Cst&&fa>=0
            fx=(-cos(C))*(Ktc*fa+Kte)*Dl+(-sin(C))*(Krc*fa+Kre)*Dl;
            fy=( sin(C))*(Ktc*fa+Kte)*Dl+(-cos(C))*(Krc*fa+Kre)*Dl;
        else
            fx=0;fy=0;
        end
        Fx=Fx+fx;Fy=Fy+fy;%%�ۼ�ÿһ����Ƭ��ÿ�������ܵ�������
        end
        FX(NB-m+1,i)=Fx;FY(NB-m+1,i)=Fy;
        Fx=0;Fy=0;%%���������ۼӱ�������
    end
    F(1,i)=sum(FX(:,i));F(2,i)=sum(FY(:,i));%%�ھ����д洢������
    %%xy���������Ҫͬʱ����
    if  i==1 %%�ж��Ƿ��ǳ�ʼ�㣨��һ���㣩
        Q1=zeros(2*NB,1);
        Q2=zeros(2*NB,1);
    else if i==2%%���ǳ�ʼ��
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
            fxs=(-cos(C))*(Ktc*fa+Kte)*Dl+(-sin(C))*(Krc*fa+Kre)*Dl;
            fys=( sin(C))*(Ktc*fa+Kte)*Dl+(-cos(C))*(Krc*fa+Kre)*Dl;
            else
                fxs=0;fys=0;
            end
            Fxs=Fxs+fxs;Fys=Fys+fys;
        end
        FXs(NB-m+1,i)=Fxs;FYs(NB-m+1,i)=Fys;
        Fxs=0;Fys=0;%%���������ۼӱ�������
    end
    Fs(1,i)=sum(FXs(:,i));Fs(2,i)=sum(FYs(:,i));%%�ھ����д洢������
end
%%��������ͼ����ɫX�򣬺�ɫY����ɫZ��
% figure(2)
% plot(Dt:Dt:Cn*Ns*Dt,F(1,:),'bo');
% hold on;
% plot(Dt:Dt:Cn*Ns*Dt,F(2,:),'ro');
% hold on;

%%���ĶԱȼӽ���%%
n1=5*Ns;
n2=9*Ns;
n3=4*Ns-200;
n4=8*Ns-200;
figure(3)
% title('Measured and predicted cutting forces','Fontsize',20)
subplot(2,1,1)
plot(Fm(n3:n4,1),FR(1,n1:n2),'m-','Markersize',7,'Markerface','white','linewidth',3.0);%%����Ҫע��F��3��N,��ȡ��������N��4
hold on;
plot(Fm(n3:n4,1),Fd(1,n1:n2),'r-','Markersize',7,'Markerface','white','linewidth',1.0);%%����Ҫע��F��3��N,��ȡ��������N��4
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
plot(Fm(n3:n4,1),FR(2,n1:n2),'m-','Markersize',7,'Markerface','white','linewidth',3.0);%%����Ҫע��F��3��N,��ȡ��������N��4
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
plot(Fm(n3:n4,1),FR(1,n1:n2),'m-','Markersize',7,'Markerface','white','linewidth',3.0);%%����Ҫע��F��3��N,��ȡ��������N��4
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
plot(Fm(n3:n4,1),FR(2,n1:n2),'m-','Markersize',7,'Markerface','white','linewidth',3.0);%%����Ҫע��F��3��N,��ȡ��������N��4
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

