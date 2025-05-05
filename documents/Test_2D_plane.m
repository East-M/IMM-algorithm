clear all;
clc;
clf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 测试两座标雷达对一个直线飞行的飞机类目标的跟踪滤波
% 雷达参数：两座标圆周扫描雷达,理论上雷达应测距测角,这里给一个简化版例子,认为x,y方向的测量是独立的
% 目标参数：匀速或匀加速飞行，初始距离\初始速度\初始加速均可设
% 滤波模型：Kalman滤波（匀速模型，即CV模型）
% 注意事项：注意过程噪声的设置，设置太小，以预测为主，比较平滑，但也有可能发散；设置太大则以量测为主，平滑效果不好
% 作业要求:
% (1) 理解代码,运行观察效果
% (2) 改变过程噪声sig_w的取值,调大或调小,观察滤波的平滑程度,并解释原因
% (3) 改变初时加速度ax或ay的取值,例如ax=0.1,观察滤波器是否发散,并解释原因
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 目标和雷达参数
T = 1;     % 扫描周期s
k = 1000;  % 观测点个数
x0 = 0;    % 目标初始位置(0km,500km)
y0 = 500e3;
vx = 300;  % 目标初时速度(m\s)
vy = 1;
ax = 0;    % 目标初时加速度(m\s^2)
ay = 0;
sig_x = 100;            % 测距精度m
sig_y = 100;            % 测角精度rad\
t = (0:k-1) * T;        % 时间集合

% 产生真实轨迹
x=  x0 + vx*t + 0.5*ax*t.^2;
y = y0 + vy*t + 0.5*ay*t.^2;

% 生成雷达原始测量数据(measure)
xm = x  + normrnd(0,sig_x,1,k);
ym = y  + normrnd(0,sig_y,1,k);

% Kalman滤波参数
n = 4;              % 状态变量维数(x xv y yv)
m = 2;              % 量测变量维数(xm,ym)
Phi = [ 1 T 0 0;    % 转移矩阵4x4,也叫做预测矩阵
        0 1 0 0;
        0 0 1 T; 
        0 0 0 1 ];    
G = [ T/2 0;        % 噪声矩阵
      1   0;
      0 T/2; 
      0    1]; 
H=[ 1 0 0 0;        % 量测矩阵,[mxn]维
    0 0 1 0];  
sig_w = [1e-2];     % 过程噪声控制参数,该参数非常重要,调节它观察效果.
Q = [sig_w^2,    0; % 过程噪声协方差
    0,   sig_w^2];
R = [sig_x^2        0; 
     0         sig_y^2]; %量测协方差
X0 = [x(2)  (x(2)-x(1))/T   y(2) (y(2)-y(1))/T]'; %滤波初始值,两点差分得到速度
P0 = [ sig_x^2        sig_x^2/T    0    0;        % 初始滤波协方差,根据差分公式得到,认为xy方向独立分布
      sig_x^2/T     2*sig_x^2/T^2  0    0; 
      0    0    sig_y^2      sig_y^2/T; 
      0    0    sig_y^2/T    2*sig_y^2/T^2];
Z = [xm; ym]; %量测矢量

%% Kalman滤波核心计算模块
for i=1:k
   [X1,P1] = Fun_KF_Predict(X0,P0,Phi,Q,G);     % 预测
   [XX,PP] = Fun_KF_Update(X1,P1,Z(:,i),H,R);   % 更新(真实的雷达数据处理需要在这两个函数中间判断点迹z是否落入预测波门,即进行数据关联,这里为简化起见,不考虑杂波和虚警)
    
    X0=XX; %为了保证循环,需要重新附值
    P0=PP;
    
    %储存结果,各个方向的滤波值
    xf(i)  = XX(1,1); 
    xvf(i) = XX(2,1);
    yf(i)  = XX(3,1);
    yvf(i) = XX(4,1);
    
    px(i) = PP(1,1);
    py(i) = PP(3,3);
end

%% 以下开始输出结果
figure(1)
% 画场景图
plot(x/1e3,y/1e3,':','LineWidth',1);   % 真实
hold on;
plot(xm/1e3,ym/1e3,'g.','LineWidth',2);% 量测
plot(xf/1e3,yf/1e3,'k','LineWidth',2); % 滤波
hold off;
xlabel('X(km)');
ylabel('Y(km)');
title('滤波航迹图');
legend('真实','观测','滤波');
axis tight

figure(2)
% 雷达P显显示
rou =sqrt( x.^2+y.^2);
theta = atan2(y,x);
rou_m = sqrt(xm.^2+ym.^2);
theta_m = atan2(ym,xm);
rou_f = sqrt(xf.^2+yf.^2);
theta_f = atan2(yf,xf);
polar(theta,rou,':');%真实
hold on
polar(theta_m,rou_m,'g.');%真实
polar(theta_f,rou_f,'k');%真实
hold off
title('雷达PPI显示');

figure(3)
% 画x方向滤波误差
plot(t,xm-x,'g:','LineWidth',2);
hold on
plot(t,xf-x,'k','LineWidth',2);
plot(t,3*sqrt(px),'r','LineWidth',2);
plot(t,-3*sqrt(px),'r','LineWidth',2);
hold off
legend('观测误差','滤波误差','3\sigma估计误差限');
xlabel('时间(s)');
ylabel('X位置误差(m)');
title('X方向滤波误差');

figure(4)
% 画x方向估计的速度
plot(t,xvf,'k-','linewidth',2); %估机值
hold on
plot(t,vx,'b'); %真实值
hold off
xlabel('时间(s)');
ylabel('速度(m/s)');
title('X方向估计速度(m/s)');

figure(5)
% 画y方向估计的速度
plot(t,yvf,'k-','linewidth',2); %估机值
hold on
plot(t,vy,'b'); %真实值
hold off
xlabel('时间(s)');
ylabel('速度(m/s)');
title('Y方向估计速度(m/s)');

