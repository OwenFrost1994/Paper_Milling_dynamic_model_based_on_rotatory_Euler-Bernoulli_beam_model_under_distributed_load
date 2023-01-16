clear;
clc;
close all;
%%ר�ż�������ģ̬����
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

DL=DL*1000;
L0=L0*1000;
%%mesh�������Ʊ�ƽ��
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
xlabel('�� order')
ylabel('L(mm)')
zlabel('Displacement(mm)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
set(get(gca,'ZLabel'),'Fontsize',20)
legend(strcat('��_0=',num2str(Wn(1,1)),'Hz'),strcat('��_1=',num2str(Wn(2,2)),'Hz'),strcat('��_2=',num2str(Wn(3,3)),'Hz'),strcat('��_3=',num2str(Wn(4,4)),'Hz'),strcat('��_4=',num2str(Wn(5,5)),'Hz'),strcat('��_5=',num2str(Wn(6,6)),'Hz'))
view(27.5, 60);

SS=1000;%%����ת��rpm
w=2*pi*SS/60;%%���߽��ٶ�
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
xlabel('�� order')
ylabel('L(mm)')
zlabel('Displacement(mm)')
set(gca, 'FontSize', 20)
set(get(gca,'XLabel'),'Fontsize',20)
set(get(gca,'YLabel'),'Fontsize',20)
set(get(gca,'ZLabel'),'Fontsize',20)
legend(strcat('��_0=',num2str(Wnw(1,1)),'Hz'),strcat('��_1=',num2str(Wnw(2,2)),'Hz'),strcat('��_2=',num2str(Wnw(3,3)),'Hz'),strcat('��_3=',num2str(Wnw(4,4)),'Hz'),strcat('��_4=',num2str(Wnw(5,5)),'Hz'),strcat('��_5=',num2str(Wnw(6,6)),'Hz'))
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
legend(strcat('��_r0=',num2str(Wnw(1,1)),'Hz'),strcat('��_s0=',num2str(Wn(1,1)),'Hz'))

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
legend(strcat('��_r1=',num2str(Wnw(2,2)),'Hz'),strcat('��_s1=',num2str(Wn(2,2)),'Hz'))

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
legend(strcat('��_r2=',num2str(Wnw(3,3)),'Hz'),strcat('��_s2=',num2str(Wn(3,3)),'Hz'))

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
legend(strcat('��_r3=',num2str(Wnw(4,4)),'Hz'),strcat('��_s3=',num2str(Wn(4,4)),'Hz'))

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
legend(strcat('��_r4=',num2str(Wnw(5,5)),'Hz'),strcat('��_s4=',num2str(Wn(5,5)),'Hz'))

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
legend(strcat('��_r5=',num2str(Wnw(6,6)),'Hz'),strcat('��_s5=',num2str(Wn(6,6)),'Hz'))