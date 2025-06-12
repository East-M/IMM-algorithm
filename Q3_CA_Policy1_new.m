%% ====================================================================== %
% Sun Yat-sen University 电子与通信工程学院
% 统计信号处理课程1班 - 第3小组
% @author: 成先锋 莫晓东 陈立邦 成泽宇
% @email:  moxd3@mail2.sysu.edu.cn
% @GitHub: https://github.com/East-M
% date:2025年5月5日
%
% #code: 大机动飞机目标的跟踪滤波技术性能仿真
% 解决问题三：基于CA模型，结合曲率挠率对数据进行滤波，并做可视化处理。
% 采用策略一：根据曲率kappa和扰率tau  动态调整Q_ca。
% 调整的参数十分重要，我们先使用真实值速度计算加速度，加加速度，进而计算得到曲率和扰率
% 然后优先调曲率的权重因子，调好后再调扰率的权重因子
% 得到合适的因子后，再带入到实时计算中，调整由实时滤波结果得到的曲率扰率。
%% ====================================================================== %
clear;
clc;
close all;

addpath utils\ data\
%% NOTICE
disp("【大机动飞机目标的跟踪滤波技术性能仿真】");
disp("解决问题一的利用CA模型对数据进行滤波，并做可视化处理。");
disp(" ");
%% 路径定义
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

% 计算加速度和加加速度，用于计算曲率和挠率
dt = T;                     % 时间间隔

%% 系统模型
% --------------------------------------------- %
%     状态方程：x(k+1) = F * x(k) + G * W(k)    |
%     量测方程：  z(k) = H * x(k) + V(k)        |
% --------------------------------------------- %

%% 初始值
X0 = [X_t(1); (X_t(2)-X_t(1))/T; (X_t(3)-2*X_t(2)+X_t(1))/T^2; % 初始位置、速度估计和加速度估计
          Y_t(1); (Y_t(2)-Y_t(1))/T; (Y_t(3)-2*Y_t(2)+Y_t(1))/T^2;
          Z_t(1); (Z_t(2)-Z_t(1))/T; (Z_t(3)-2*Z_t(2)+Z_t(1))/T^2];
P0 = diag([sigma_r^2, 50^2, 25^2, sigma_r^2, 50^2, 25^2, sigma_r^2, 50^2, 25^2]);
xNum = 9;          % 状态矢量参数量
zNum = 9;          % 观测矢量参数量

%% CA模型参数
% 状态转移矩阵
F_ca = blkdiag([1,T,T^2/2;0,1,T;0,0,1],[1,T,T^2/2;0,1,T;0,0,1],[1,T,T^2/2;0,1,T;0,0,1]);

% 系统扰动矩阵
G_ca = blkdiag([T^3/6;T^2/2;T],[T^3/6;T^2/2;T],[T^3/6;T^2/2;T]);

% 量测矩阵
% H_ca = blkdiag([1,0,0],[1,0,0],[1,0,0]);
H_ca = eye(9);

% 过程噪声协方差矩阵
Q_ca = sigma^2 * (G_ca * G_ca');

% 量测噪声协方差矩阵
R_ca = sigma_r^2 * eye(zNum);

%% 蒙特卡洛仿真实验
X_ca_filted = zeros(xNum,N);
rmse_obs = 0;
cost_time = zeros(Num,1);

disp("[NOTICE] 正在进行蒙特卡洛实验...");
for num = 1:Num
    % 叠加量测误差，生成观测数据
    x = X_t + normrnd(0, sigma_r, size(X_t));
    y = Y_t + normrnd(0, sigma_r, size(Y_t));
    z = Z_t + normrnd(0, sigma_r, size(Z_t));
    vx = XV_t + normrnd(0, sigma_r, size(XV_t));
    vy = YV_t + normrnd(0, sigma_r, size(YV_t));
    vz = ZV_t + normrnd(0, sigma_r, size(ZV_t));
    ax = [0; diff(XV_t)] + normrnd(0, sigma_r, size(XV_t));
    ay = [0; diff(YV_t)] + normrnd(0, sigma_r, size(YV_t));
    az = [0; diff(ZV_t)] + normrnd(0, sigma_r, size(ZV_t));
    
    rmse_obs = rmse_obs + sqrt(mean((x(obs_idx_beg:end)-X_t(obs_idx_beg:end)).^2 + ...
                         (y(obs_idx_beg:end)-Y_t(obs_idx_beg:end)).^2 + ...
                         (z(obs_idx_beg:end)-Z_t(obs_idx_beg:end)).^2));

    % 初始化变量
    Z = [x';vx';ax';y';vy';ay';z';vz';az'];       % 观测数据
    Pk = P0;
    I = eye(size(P0));
    kappa = 0;
    tau = 0;
    alpha = 1000000;
    beta = 1000;
    X_pre = X0;

    tic;

    for k = 1:N
        % 在线计算卡尔曼增益 以及 预测&滤波的误差协方差矩阵
        
        % 先计算kappa 和 tau
        
        if k >= 3
            v = [X_ca_filted(2,k-1), X_ca_filted(5,k-1), X_ca_filted(8,k-1)]; % 速度向量
            a = [(X_ca_filted(2,k-1) - X_ca_filted(2,k-2)) / dt, 
                (X_ca_filted(5,k-1) - X_ca_filted(5,k-2)) / dt,
                (X_ca_filted(8,k-1) - X_ca_filted(8,k-2)) / dt];  % 加速度向量
            j = [(X_ca_filted(3,k-1) - X_ca_filted(3,k-2)) / dt, 
                (X_ca_filted(6,k-1) - X_ca_filted(6,k-2)) / dt,
                (X_ca_filted(9,k-1) - X_ca_filted(9,k-2)) / dt];  % 加加速度向量
        else
            v = [X_pre(2), X_pre(5), X_pre(8)];     % 速度向量
            a = [X_pre(3), X_pre(6), X_pre(9)];     % 加速度向量
            j = [0,0,0];
        end
        [kappa, tau] = kappa_tau_cal(v,a,j);        % 计算得到此时的曲率和扰率，由于是cv模型，加速度始终是0，所以tau始终是0

        if kappa ~= 0 || tau ~=0
            scale_Q = 1 + alpha * kappa + beta * tau;% alpha beta 为调节系数
        else
            scale_Q = 1;
        end
        Q_ca_adaptive = scale_Q * Q_ca;             % 调节过程噪声矩阵

        Pk_estm = F_ca * Pk * F_ca' + Q_ca_adaptive;                % 预测
        K_ca = Pk_estm * H_ca' / (H_ca * Pk_estm * H_ca' + R_ca);   % 卡尔曼增益
        Pk = (I - K_ca * H_ca) * Pk_estm;                           % 更新

        X_pre = KalmanFilter(Z(:,k), X_pre, F_ca, H_ca, K_ca);
        X_ca_filted(:,k) = X_pre + X_ca_filted(:,k);
    end
    cost_time(num) = toc;
    % disp("第"+num2str(num)+"次蒙特卡洛实验结束.");
end
disp("【蒙特卡洛实验结束!】 共进行"+num2str(Num)+"次实验。");
X_ca_filted = X_ca_filted ./ Num;
rmse_obs = rmse_obs ./ Num;
cost_time = mean(cost_time);
disp("平均实验耗时：" + num2str(cost_time*1000) + " ms");

x_filted = X_ca_filted(1,:);
xv_filted = X_ca_filted(2,:);
y_filted = X_ca_filted(4,:);
yv_filted = X_ca_filted(5,:);
z_filted = X_ca_filted(7,:);
zv_filted = X_ca_filted(8,:);

% 位置估计RMSE
rmse_xyz_ca = sqrt(mean((x_filted(obs_idx_beg:end)' - X_t(obs_idx_beg:end)).^2 + ...
                        (y_filted(obs_idx_beg:end)' - Y_t(obs_idx_beg:end)).^2 + ...
                        (z_filted(obs_idx_beg:end)' - Z_t(obs_idx_beg:end)).^2));
disp("=== 观测数据平均位置RMSE："+num2str(rmse_obs)+" m");
disp("=== 卡尔曼滤波后位置RMSE："+num2str(rmse_xyz_ca)+" m");

%% CA模型可视化滤波结果
vis = true;
if vis == true
    figure('name',"实验数据",'Position', [100, 100, 1200, 800]);
    plot3(X_t, Y_t, Z_t, 'LineWidth', 2); title("飞机真实轨迹");
    xlabel("x方向 (m)"); ylabel("y方向 ((m)"); zlabel("z方向 (m)");
    grid on; hold on;
    plot3(x, y, z, '.', 'MarkerSize', 8); title("飞行轨迹观测数据");
    % plot3(x_filted, y_filted, z_filted, 'LineWidth', 2);
    legend("Real Trajectory", "Observed Data");
    
    figure('name',"滤波结果3维图",'Position', [100, 100, 1200, 800]);
    plot3(x, y, z, '.', 'MarkerSize', 8);
    title("CA with kappa and tau");
    xlabel("x(m)"); ylabel("y(m)"); zlabel("z(m)");
    grid on; hold on;
    plot3(x_filted, y_filted, z_filted, 'LineWidth', 2);
    plot3(X_t, Y_t, Z_t, 'green', 'LineWidth', 2);
    legend("Observed Data", "Result", "Real Trajectory");
    exportgraphics(gcf, "策略一：CA+kappa+tau滤波结果.pdf", "ContentType", "vector")
    

    figure('name',"CA模型的速度卡尔曼滤波结果",'Position', [1300, 100, 1200, 800]);
    subplot(3,1,1);
    plot(XV_t); xlabel("k"); ylabel("m/s"); title("xv");
    hold on; grid on;
    plot(xv_filted);
    legend("True Speed","Result");
    
    subplot(3,1,2);
    plot(YV_t); xlabel("k"); ylabel("m/s"); title("yv");
    hold on; grid on;
    plot(yv_filted);
    legend("True Speed","Result");
    
    subplot(3,1,3);
    plot(ZV_t); xlabel("k"); ylabel("m/s"); title("zv");
    hold on; grid on;
    plot(zv_filted);
    legend("True Speed","Result");
    exportgraphics(gcf, "策略一：CA+kappa+tau模型的速度的滤波结果.pdf", "ContentType", "vector")

    
    figure('name',"CA模型的位置卡尔曼滤波结果",'Position', [1300, 100, 1200, 800]);
    subplot(3,1,1);
    plot(X_t); xlabel("k"); ylabel("m"); title("x");
    hold on; grid on;
    plot(x,'.');
    plot(x_filted);
    legend("Real Trajectory","Observed Data","Result");
    
    subplot(3,1,2);
    plot(Y_t); xlabel("k"); ylabel("m"); title("y");
    hold on; grid on;
    plot(y,'.');
    plot(y_filted);
    legend("Real Trajectory","Observed Data","Result");
    
    subplot(3,1,3);
    plot(Z_t); xlabel("k"); ylabel("m"); title("z");
    hold on; grid on;
    plot(z,'.');
    plot(z_filted);
    legend("Real Trajectory","Observed Data","Result");
    exportgraphics(gcf, "策略一：CA+kappa+tau模型的位置的滤波结果.pdf", "ContentType", "vector")



    disp("可视化结果完成.");
end