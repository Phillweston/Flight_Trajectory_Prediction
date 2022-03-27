%---------------------------------------------%
%					      %
%          工作室提供代做matlab仿真	      %
%					      %
%  详情请访问：http://cn.mikecrm.com/5k6v1DP  %
%					      %
%---------------------------------------------%

%使用卡尔曼滤波方法对飞行航班进行轨迹预测
%数据来源：FlightAware（https://zh.flightaware.com）
%航    班：CES9937     宁波栎社国际机场飞往成都双流国际机场
%飞行时间：2018-07-17  07:07-09:55
%说明：取起飞后,前20组数据作为实验数据。对时间点进行近似取值，假设每隔30s上报一次数据
clear;
clc;
%采样点的个数
N=228;
%测试数据：纬度
latitude=load('C:\workspace\matlab workspace\CES9937\latitude.txt');
%真实维度值
lat=latitude;
%卡尔曼滤波处理的状态，即估计值
lat_kf=zeros(1,N);
%测报值
lat_z=zeros(1,N);
P=zeros(1,N);
%初始纬度值
lat(1)=29.8131;
%初始值的协方差
P(1)=0.09;
%初始测报值
lat_z(1)=29.8027;
%初始估计状态。假设和初始测报值相同
lat_kf(1)=lat_z(1);
%噪声方差
%系统噪声方差
Q=0.1;
%测量噪声方差
R=0.001;
%方差决定噪声大小
W=sqrt(Q)*randn(1,N);
V=sqrt(R)*randn(1,N);
%系统矩阵
F=1;
G=1;
H=1;
%本系统状态为1维
I=eye(1);
%模拟纬度测报，并滤波
for k=2:N
    %随时间推移，飞行纬度逐渐变化
    %k时刻的真是纬度值是测报仪器不知道的，测报值可能是无限接近于真实值，但并不是真实值
    %lat(k)=F*lat(k-1)+G*W(k-1);
    %纬度在k时刻的测报值
    lat_z(k)=H*lat(k)+V(k);
    %kalman滤波
    %有了k时刻的测报值lat_z(k)和k-1时刻的状态，那么就可以进行滤波了
    %状态预测
    lat_pre=F*lat_kf(k-1);
    %协方差预测
    P_pre=F*P(k-1)*F'+Q;
    %计算卡尔曼增益
    Kg=P_pre*inv(H*P_pre*H'+R);
    %新息
    e=lat_z(k)-H*lat_pre;
    %状态更新
    lat_kf(k)=lat_pre+Kg*e;
    %协方差更新
    P(k)=(I-Kg*H)*P_pre;
end
%计算误差
%测量值与真实值之间的偏差
Err_Messure=zeros(1,N);
%kalman估计与真实值之间的偏差
Err_Kalman=zeros(1,N);
for k=1:N
    Err_Messure(k)=abs(lat_z(k)-lat(k));
    Err_Kalman(k)=abs(lat_kf(k)-lat(k));
end
t=1:N;
%滤波效果图
figure
plot(t,lat,'g-',t,lat_z,'b-',t,lat_kf,'r-');
legend('真实值','观测值','kalman滤波值');
xlabel('测报时间点');
ylabel('纬度值');
%误差分析图
figure
plot(t,Err_Messure,'b-',t,Err_Kalman,'r-');
legend('测报偏差','kalman滤波偏差');
xlabel('测报时间点');
ylabel('纬度偏差值');