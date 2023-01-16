clear;
clc;
close all;
%%������ѧϵͳ���������µķ��棬���ߵ�ЧΪ���ξ�̬��ģ�ͣ������%%
%%����������ģ�͵Ļ����������������%%
DL=0.02;%%����Ԫ�ĳ���,mm
E=719000;%%����ģ��MPa=N/mm^2
p=14400e-6;%%�ܶ�kg/(cm^3)
L0=70;%%70����������mm��0�Žڵ㲻�ùܣ�ֵ��֪
L1=35;%%��һ�ν���ĳ��ȣ�mm
L2=L0-L1;%%�ڶ��ν���ĳ��ȣ�mm
D1=12.7;%%��һ��ֱ����mm
D2=10;%%�ڶ���ֱ����mm
I1=pi*(D1)^4/64;%%��һ�ι��Ծ�,mm^4
S1=pi*(D1)^2/4;%%��һ�ν����,mm^2
I2=pi*(D2)^4/64;%%�ڶ��ι��Ծ�,mm^4
S2=pi*(D2)^2/4;%%�ڶ��ν����,mm^2
NB=ceil(L0/DL);%%����Ԫ�ĸ���

%%���Ծغͽ������ɢ%%
for i=1:1:NB
    l(i)=i*DL-0.5*DL;%%�����Ƭ���е�λ��
    if l(i)<=L1
        I(i)=I1;
        S(i)=S1;
    else
        I(i)=I2;
        S(i)=S2;
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
%%�������ÿ��ʱ������ϵ�����󣬳����������ջ�þ�̬�նȾ������������%%
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
    K(i,:)=E*(ID(i,1)*A1(i,:)/(12*DL)+ID(i,2)*A2(i,:)/(2*DL^2)+ID(i,3)*A3(i,:)/(DL^3));
end
%%�����������%%
for i=1:1:NB
    M(i,i)=p*DL*S(i)*1e-6;
end

%%��μ���������
%%�����ǵ���ƫ�ĺ���б%%

%%���߼��β���%%
D=D1;%%���߰뾶
N=2;%%���߳���
B=pi/6;%%����������
Cp=2*pi/N;%%�ݼ��
%%΢Ԫ���ȼ�Ϊ��������Ԫ��΢Ԫ����

%%���ϲ���%%����AL7075%%
% Ktc=951.751;%%���������ϵ��
% Kte=14.0371;%%�����п���ϵ��
% Krc=608.561;%%���������ϵ��
% Kre=16.5002;%%�����п���ϵ��
% Kac=288.478;%%���������ϵ��
% Kae=-1.25118;%%�����п���ϵ��
%%���ϲ���%%����AL6061-T6511%%
Ktc=613.92;%%���������ϵ��
Kte=12.785;%%�����п���ϵ��
Krc=330.75;%%���������ϵ��
Kre=12.902;%%�����п���ϵ��
Kac=78.03;%%���������ϵ��
Kae=1.3672;%%�����п���ϵ��

%%�ӹ�����%%
%%ϳ����ʽ:˳ϳ%%
Cm=1;%%ϳ����ʽ��˳ϳΪ1����ϳΪ0
S=2000;%%����ת��rmp
f=400;%%�����ٶ�mm/min
fs=10000;%%����Ƶ��
ap=2;%%��������mm
ae=D;%%��������mm
Cn=5;%%Ȧ��circle number

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
Dl=DL;
FX=zeros(NB,Ns*Cn);
FY=zeros(NB,Ns*Cn);
FZ=zeros(NB,Ns*Cn);%%��¼ÿһ������Ԫ��ÿһʱ�̵���������ʵ����ֻ��ͷ��������Ԫ����
Fx=0;
Fy=0;
Fz=0;
F=zeros(3,Ns*Cn);%%�洢���������������
x=zeros(NB,Ns*Cn);
y=zeros(NB,Ns*Cn);
Dx=zeros(ap/Dl,Ns*Cn);%%������x�����λ�Ʊ仯����fv�ĵ�һ������,���ﶯ̬��м���ֻ���������εĲ��֣���ŷ�ʽ���б�ŷ�ʽһ��
Dy=zeros(ap/Dl,Ns*Cn);%%������y�����λ�Ʊ仯����fv�ĵڶ�������
apx=zeros(1,Cn*N);
apy=zeros(1,Cn*N);
K1=(M/(Dt^2)+K);%%΢�ַ�����ɢ�����
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
            Dx(m,i)=(x(NB-m+1,i)-0);%%ע�⵽����ѧ�����������е�λΪmm
            Dy(m,i)=(y(NB-m+1,i)-0);%%ע��ʱ����ÿһ�У���ͬ���б������
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
            Dx(m,i)=(x(i)-min(apx(1:q)));
            Dy(m,i)=(y(i)-min(apy(1:q)));
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
            fz=(Kac*fa+Kae)*Dl;
        else
            fx=0;fy=0;fz=0;
        end
        Fx=Fx+fx;Fy=Fy+fy;Fz=Fz+fz;%%�ۼ�ÿһ����Ƭ��ÿ�������ܵ�������
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
plot(Dt:Dt:Cn*Ns*Dt,F(1,:),'b-','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
xlabel('time(s)')
ylabel('Force(N)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
plot(Dt:Dt:Cn*Ns*Dt,F(2,:),'r-','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
xlabel('time(s)')
ylabel('Force(N)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
plot(Dt:Dt:Cn*Ns*Dt,F(3,:),'g-','Markersize',7,'Markerface','white','linewidth',3.0);
grid on;
xlabel('time(s)')
ylabel('Force(N)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
title('Summary of forces imposed on cutter');
legend('X predicted','Y predicted','Z predicted')

figure(6)
surf(Dt:Dt:Dt*Ns*Cn,DL:DL:L0,x)
grid on;
hold on;
xlabel('time(s)')
ylabel('Length(mm)')
zlabel('Displacement in x direction(mm)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
set(get(gca,'ZLabel'),'Fontsize',20)
title('Vibration displacement in x');
figure(7)
surf(Dt:Dt:Dt*Ns*Cn,DL:DL:L0,y)
grid on;
hold on;
xlabel('time(s)')
ylabel('Length(mm)')
zlabel('Displacement in y direction(mm)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
set(get(gca,'ZLabel'),'Fontsize',20)
title('Vibration displacement in y');

% figure(8)
% surf(Dt:Dt:Dt*Ns*Cn,Dl:Dl:ap,Dx)
% grid on;
% hold on;
% xlabel('time(s)')
% ylabel('Length(mm)')
% zlabel('Dynamic displacement difference in x direction(mm)')
% set(gca, 'FontSize', 20)
% set(get(gca,'XLabel'),'Fontsize',20)
% set(get(gca,'YLabel'),'Fontsize',20)
% set(get(gca,'ZLabel'),'Fontsize',20)
% title('Dynamic displacement difference in x');
% 
% figure(9)
% surf(Dt:Dt:Dt*Ns*Cn,Dl:Dl:ap,Dy)
% grid on;
% hold on;
% xlabel('time(s)')
% ylabel('Length(m)')
% zlabel('Dynamic displacement difference in y direction(mm)')
% set(gca, 'FontSize', 20)
% set(get(gca,'XLabel'),'Fontsize',20)
% set(get(gca,'YLabel'),'Fontsize',20)
% set(get(gca,'ZLabel'),'Fontsize',20)
% title('Dynamic displacement difference in y');

figure(10)
subplot(3,1,1)
plot(Dt:Dt:Dt*Ns*Cn,x(NB,:),'r+-','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
plot(Dt:Dt:Dt*Ns*Cn,y(NB,:),'b+-','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
grid on;
xlabel('time(s)')
ylabel('Displacement(mm)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
set(get(gca,'ZLabel'),'Fontsize',20)
title('Dynamic displacement of note 3500');
legend('Displacement in x','Displacement in y')
subplot(3,1,2)
plot(Dt:Dt:Dt*Ns*Cn,x(NB-250,:),'r+-','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
plot(Dt:Dt:Dt*Ns*Cn,y(NB-250,:),'b+-','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
grid on;
xlabel('time(s)')
ylabel('Displacement(mm)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
set(get(gca,'ZLabel'),'Fontsize',20)
title('Dynamic displacement of note 3250');
legend('Displacement in x','Displacement in y')
subplot(3,1,3)
plot(Dt:Dt:Dt*Ns*Cn,x(NB-500,:),'r+-','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
plot(Dt:Dt:Dt*Ns*Cn,y(NB-500,:),'b+-','Markersize',7,'Markerface','white','linewidth',3.0);
hold on;
grid on;
xlabel('time(s)')
ylabel('Displacement(mm)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
set(get(gca,'ZLabel'),'Fontsize',20)
title('Dynamic displacement of note 3000');
legend('Displacement in x','Displacement in y')