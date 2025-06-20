%% ====================================================================== %
% Sun Yat-sen University 电子与通信工程学院
% 统计信号处理课程1班 - 第3小组
% @author: 成先锋 莫晓东 陈立邦 成泽宇
% @email:  moxd3@mail2.sysu.edu.cn
% @GitHub: https://github.com/East-M
% date:2025年5月5日
%
% #code: 大机动飞机目标的跟踪滤波技术性能仿真
% 解决问题一的利用Singer模型对数据进行滤波，并做可视化处理。
%% ====================================================================== %
clear;
clc;
close all;

addpath utils\ data\
%% NOTICE
disp("【大机动飞机目标的跟踪滤波技术性能仿真】");
disp("解决问题一的利用Singer模型对数据进行滤波，并做可视化处理。");
disp(" ");
%% 路径定义
K_singer_Path = "./data/K_singer.mat";
DATA_path = "./data/dataENU.txt";

%% 读取数据
format = '%f %f %f %f %f %f %f';
fid = fopen(DATA_path, 'r');
data = textscan(fid, format);
data{1} = linspace(0, 194, 389)';
fclose(fid);

%% 参数设置
Num = 100;       % 蒙特卡洛仿真次数
T = 0.5;         % 扫描周期
N = 389;         % 观测点个数
sigma = 10;           % 过程噪声方差(为设计参数，需要反复调试)
sigma_r = 100;        % 量测误差，单位: m
obs_idx_beg = 1;     % 计算RMSE的起始位置（因为卡尔曼滤波需要经过一段时间收敛后，滤波效果才好，前面滤波得到的结果不是很准确，会影响RMSE的计算）

%% 真实轨迹数据
X_t = data{2};     % 目标位置
Y_t = data{3};
Z_t = data{4};
XV_t = data{5};    % 目标速度
YV_t = data{6};
ZV_t = data{7};

%% 系统模型
% --------------------------------------------- %
%     状态方程：x(k+1) = F * x(k) + G * W(k)    |
%     量测方程：  z(k) = H * x(k) + V(k)        |
% --------------------------------------------- %
%% 初始值
xNum = 9;          % 状态矢量参数量
zNum = 3;          % 观测矢量参数量
P0 = diag([sigma_r^2, 50^2, 25^2, sigma_r^2, 50^2, 25^2, sigma_r^2, 50^2, 25^2]);
X0 = [X_t(1); (X_t(2)-X_t(1))/T; (X_t(3)-2*X_t(2)+X_t(1))/T^2; % 初始位置、速度估计和加速度估计
          Y_t(1); (Y_t(2)-Y_t(1))/T; (Y_t(3)-2*Y_t(2)+Y_t(1))/T^2;
          Z_t(1); (Z_t(2)-Z_t(1))/T; (Z_t(3)-2*Z_t(2)+Z_t(1))/T^2];

%% Singer模型参数
sigma_m = 20;       % Singer模型的加速度方差
alpha = 1/5;       % 机动频率，经验值。

% 状态转移矩阵
e = exp(1);
blk = [1,T,(alpha*T-1+e^(-alpha*T))/alpha^2;0,1,(1-e^(-alpha*T))/alpha;0,0,e^(-alpha*T)];
F_singer = blkdiag(blk,blk,blk);

% 量测矩阵
H_singer = blkdiag([1,0,0],[1,0,0],[1,0,0]);

% 过程噪声协方差矩阵
q11 = (1-e^(-2*alpha*T)+2*alpha*T+2*alpha^3*T^3/3-2*alpha^2*T^2-4*alpha*T*e^(-alpha*T))/(2*alpha^5);
q12 = (e^(-2*alpha*T)+1-2*e^(-alpha*T)+2*alpha*T*e^(-alpha*T)-2*alpha*T+alpha^2*T^2)/(2*alpha^4);
q13 = (1-e^(-2*alpha*T)-2*alpha*T*e^(-alpha*T))/(2*alpha^3);
q22 = (4*e^(-alpha*T)-3-e^(-2*alpha*T)+2*alpha*T)/(2*alpha^3);
q23 = (e^(-2*alpha*T)+1-2*e^(-alpha*T))/(2*alpha^2);
q33 = (1-e^(-2*alpha*T))/(2*alpha);

Qblk = [q11,q12,q13;q12,q22,q23;q13,q23,q33];
Q_singer = 2*alpha*sigma_m^2*blkdiag(Qblk,Qblk,Qblk);

% 量测噪声协方差矩阵
R_singer = sigma_r^2 * eye(zNum);

% ----------------- 离线计算卡尔曼增益 ----------------- %
disp("[NOTICE] 正在加载卡尔曼增益...");
if exist(K_singer_Path,"file") == 0
    % 未找到离线计算的卡尔曼增益文件，进行计算和保存
    K_singer = KalmanGain(P0, F_singer, H_singer, Q_singer, R_singer);
    save(K_singer_Path, "K_singer");
    disp("已离线计算卡尔曼增益! 保存于【"+K_singer_Path+"】");
else
    load(K_singer_Path);
    disp("已加载Singer模型的卡尔曼增益!");
end
% ------------------------------------------------------ %

%% 蒙特卡洛仿真实验
X_singer_filted = zeros(xNum,N);
rmse_obs = 0;
cost_time = zeros(Num,1);

disp("[NOTICE] 正在进行蒙特卡洛实验...");
for num = 1:Num
    % 叠加量测误差，生成观测数据
    x = X_t + normrnd(0, sigma_r, size(X_t));
    y = Y_t + normrnd(0, sigma_r, size(Y_t));
    z = Z_t + normrnd(0, sigma_r, size(Z_t));

    rmse_obs = rmse_obs + sqrt(mean((x(obs_idx_beg:end)-X_t(obs_idx_beg:end)).^2 + ...
                         (y(obs_idx_beg:end)-Y_t(obs_idx_beg:end)).^2 + ...
                         (z(obs_idx_beg:end)-Z_t(obs_idx_beg:end)).^2));

    % 初始值
    Z = [x';y';z'];      % 观测数据

    X_pre = X0;
    tic;
    for k = 1:N
        X_pre = KalmanFilter(Z(:,k), X_pre, F_singer, H_singer, K_singer(:,:,k));
        X_singer_filted(:,k) = X_pre + X_singer_filted(:,k);
    end
    cost_time(num) = toc;
    % disp("第"+num2str(num)+"次蒙特卡洛实验结束.");
end
disp("【蒙特卡洛实验结束!】 共进行"+num2str(Num)+"次实验。");
X_singer_filted = X_singer_filted ./ Num;
rmse_obs = rmse_obs ./ Num;
cost_time = mean(cost_time);
disp("平均实验耗时：" + num2str(cost_time*1000) + " ms");

x_filted = X_singer_filted(1,:);
xv_filted = X_singer_filted(2,:);
y_filted = X_singer_filted(4,:);
yv_filted = X_singer_filted(5,:);
z_filted = X_singer_filted(7,:);
zv_filted = X_singer_filted(8,:);

% 位置估计RMSE
rmse_xyz_singer = sqrt(mean((x_filted(obs_idx_beg:end)' - X_t(obs_idx_beg:end)).^2 + ...
                        (y_filted(obs_idx_beg:end)' - Y_t(obs_idx_beg:end)).^2 + ...
                        (z_filted(obs_idx_beg:end)' - Z_t(obs_idx_beg:end)).^2));
disp("=== 观测数据平均位置RMSE："+num2str(rmse_obs)+" m");
disp("=== 卡尔曼滤波后位置RMSE："+num2str(rmse_xyz_singer)+" m");

%% Singer模型可视化滤波结果
figure('name',"实验数据",'Position', [100, 100, 1200, 800]);
plot3(X_t, Y_t, Z_t, 'LineWidth', 2); title("飞机真实轨迹");
xlabel("x方向 (m)"); ylabel("y方向 ((m)"); zlabel("z方向 (m)");
grid on; hold on;
plot3(x, y, z, '.', 'MarkerSize', 8); title("飞行轨迹观测数据");
% plot3(x_filted, y_filted, z_filted, 'LineWidth', 2);
legend("真实轨迹", "观测数据");

figure('name',"滤波结果3维图",'Position', [100, 100, 1200, 800]);
plot3(x, y, z, '.', 'MarkerSize', 8);
title("Singer模型的卡尔曼滤波结果");
xlabel("x方向 (m)"); ylabel("y方向 ((m)"); zlabel("z方向 (m)");
grid on; hold on;
plot3(x_filted, y_filted, z_filted, 'LineWidth', 2);
legend("观测数据", "滤波结果");

figure('name',"Singer模型的速度卡尔曼滤波结果",'Position', [1300, 100, 1200, 800]);
subplot(3,1,1);
plot(XV_t); xlabel("k"); ylabel("m/s"); title("xv分量");
hold on; grid on;
plot(xv_filted);
legend("真实速度","滤波结果");

subplot(3,1,2);
plot(YV_t); xlabel("k"); ylabel("m/s"); title("yv分量");
hold on; grid on;
plot(yv_filted);
legend("真实速度","滤波结果");

subplot(3,1,3);
plot(ZV_t); xlabel("k"); ylabel("m/s"); title("zv分量");
hold on; grid on;
plot(zv_filted);
legend("真实速度","滤波结果");

figure('name',"Singer模型的位置卡尔曼滤波结果",'Position', [1300, 100, 1200, 800]);
subplot(3,1,1);
plot(X_t); xlabel("k"); ylabel("m"); title("x分量");
hold on; grid on;
plot(x,'.');
plot(x_filted);
legend("真实轨迹","观测数据","滤波结果");

subplot(3,1,2);
plot(Y_t); xlabel("k"); ylabel("m"); title("y分量");
hold on; grid on;
plot(y,'.');
plot(y_filted);
legend("真实轨迹","观测数据","滤波结果");

subplot(3,1,3);
plot(Z_t); xlabel("k"); ylabel("m"); title("z分量");
hold on; grid on;
plot(z,'.');
plot(z_filted);
legend("真实轨迹","观测数据","滤波结果");
disp("可视化结果完成.");