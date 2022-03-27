%---------------------------------------------%
%					      %
%          �������ṩ����matlab����	      %
%					      %
%  ��������ʣ�http://cn.mikecrm.com/5k6v1DP  %
%					      %
%---------------------------------------------%

%ʹ�ÿ������˲������Է��к�����й켣Ԥ��
%������Դ��FlightAware��https://zh.flightaware.com��
%��    �ࣺCES9937     ����������ʻ��������ɶ�˫�����ʻ���
%����ʱ�䣺2018-07-17  07:07-09:55
%˵����ȡ��ɺ�,ǰ20��������Ϊʵ�����ݡ���ʱ�����н���ȡֵ������ÿ��30s�ϱ�һ������
clear;
clc;
%������ĸ���
N=228;
%�������ݣ�γ��
latitude=load('C:\workspace\matlab workspace\CES9937\latitude.txt');
%��ʵά��ֵ
lat=latitude;
%�������˲������״̬��������ֵ
lat_kf=zeros(1,N);
%�ⱨֵ
lat_z=zeros(1,N);
P=zeros(1,N);
%��ʼγ��ֵ
lat(1)=29.8131;
%��ʼֵ��Э����
P(1)=0.09;
%��ʼ�ⱨֵ
lat_z(1)=29.8027;
%��ʼ����״̬������ͳ�ʼ�ⱨֵ��ͬ
lat_kf(1)=lat_z(1);
%��������
%ϵͳ��������
Q=0.1;
%������������
R=0.001;
%�������������С
W=sqrt(Q)*randn(1,N);
V=sqrt(R)*randn(1,N);
%ϵͳ����
F=1;
G=1;
H=1;
%��ϵͳ״̬Ϊ1ά
I=eye(1);
%ģ��γ�Ȳⱨ�����˲�
for k=2:N
    %��ʱ�����ƣ�����γ���𽥱仯
    %kʱ�̵�����γ��ֵ�ǲⱨ������֪���ģ��ⱨֵ���������޽ӽ�����ʵֵ������������ʵֵ
    %lat(k)=F*lat(k-1)+G*W(k-1);
    %γ����kʱ�̵Ĳⱨֵ
    lat_z(k)=H*lat(k)+V(k);
    %kalman�˲�
    %����kʱ�̵Ĳⱨֵlat_z(k)��k-1ʱ�̵�״̬����ô�Ϳ��Խ����˲���
    %״̬Ԥ��
    lat_pre=F*lat_kf(k-1);
    %Э����Ԥ��
    P_pre=F*P(k-1)*F'+Q;
    %���㿨��������
    Kg=P_pre*inv(H*P_pre*H'+R);
    %��Ϣ
    e=lat_z(k)-H*lat_pre;
    %״̬����
    lat_kf(k)=lat_pre+Kg*e;
    %Э�������
    P(k)=(I-Kg*H)*P_pre;
end
%�������
%����ֵ����ʵֵ֮���ƫ��
Err_Messure=zeros(1,N);
%kalman��������ʵֵ֮���ƫ��
Err_Kalman=zeros(1,N);
for k=1:N
    Err_Messure(k)=abs(lat_z(k)-lat(k));
    Err_Kalman(k)=abs(lat_kf(k)-lat(k));
end
t=1:N;
%�˲�Ч��ͼ
figure
plot(t,lat,'g-',t,lat_z,'b-',t,lat_kf,'r-');
legend('��ʵֵ','�۲�ֵ','kalman�˲�ֵ');
xlabel('�ⱨʱ���');
ylabel('γ��ֵ');
%������ͼ
figure
plot(t,Err_Messure,'b-',t,Err_Kalman,'r-');
legend('�ⱨƫ��','kalman�˲�ƫ��');
xlabel('�ⱨʱ���');
ylabel('γ��ƫ��ֵ');