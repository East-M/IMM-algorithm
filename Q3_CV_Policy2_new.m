%% ====================================================================== %
% Sun Yat-sen University 电子与通信工程学院
% 统计信号处理课程1班 - 第3小组
% @author: 成先锋 莫晓东 陈立邦 成泽宇
% date:2025年5月6日
%
% #code: 大机动飞机目标的跟踪滤波技术性能仿真
% 解决问题三：基于CV模型，结合曲率挠率对数据进行滤波，并做可视化处理。
% 采用策略二：将曲率和挠率直接当成参数值，进行扩维滤波实时估计。
%% ====================================================================== %
clear;
clc;
close all;

addpath utils\ data\
%% NOTICE
disp("【大机动飞机目标的跟踪滤波技术性能仿真】");
disp("基于CV模型，结合曲率挠率对数据进行滤波，并做可视化处理。");
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
sigma_kappa = 0.001;  % 曲率过程噪声
sigma_tau = 0.0005;   % 挠率过程噪声
obs_idx_beg = 1;     % 计算RMSE的起始位置（因为卡尔曼滤波需要经过一段时间收敛后，滤波效果才好，前面滤波得到的结果不是很准确，会影响RMSE的计算）

%% 真实轨迹数据
X_t = data{2};     % 目标位置
Y_t = data{3};
Z_t = data{4};
XV_t = data{5};    % 目标速度
YV_t = data{6};
ZV_t = data{7};
XA_t = [0; diff(XV_t)]; % 目标加速度
YA_t = [0; diff(YV_t)];
ZA_t = [0; diff(ZV_t)];

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
      Z_t(1); (Z_t(2)-Z_t(1))/T; (Z_t(3)-2*Z_t(2)+Z_t(1))/T^2;
      0; 0];
P0 = diag([sigma_r^2, 50^2, 0, sigma_r^2, 50^2, 0, sigma_r^2, 50^2, 0, ...
    100^2, 100^2]); % 设置初始位置，初始速度，初始加速度，kappa，tau的估计（滤波）误差方差

xNum = 11;         % 状态矢量参数量
zNum = 11;          % 观测矢量参数量

%% CV模型参数
% 状态转移矩阵
F_cv = blkdiag([1,T,0;0,1,0;0,0,0],[1,T,0;0,1,0;0,0,0],[1,T,0;0,1,0;0,0,0]);
F_cv_augmented = blkdiag(F_cv, eye(2));             % 假设曲率和挠率随机游走

% 系统扰动矩阵：
G_cv = blkdiag([0.5*T^2;T;0],[0.5*T^2;T;0],[0.5*T^2;T;0]);
G_cv_augmented = blkdiag(G_cv, eye(2));             % 曲率和挠率随机游走下，过程噪声扩展

% 量测矩阵
% H_cv = blkdiag([1,0,0;0,1,0],[1,0,0;0,1,0],[1,0,0;0,1,0]);
% H_cv_augmented = [H_cv, zeros(6,2)];                % 仅观测位置
H_cv = blkdiag([1,0,0;0,1,0;0,0,0],[1,0,0;0,1,0;0,0,0],[1,0,0;0,1,0;0,0,0]);
H_cv_augmented = blkdiag(H_cv, eye(2));                % 仅观测位置

% 过程噪声协方差矩阵
Q_cv = sigma^2 * (G_cv * G_cv');
Q_cv_augmented = blkdiag(Q_cv, diag([sigma_kappa^2, sigma_tau^2])); % 新增噪声项

% 量测噪声协方差矩阵
R_cv = sigma_r^2 * eye(zNum);

%% 蒙特卡洛仿真实验
X_cv_filted = zeros(xNum,N);
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
    ax = zeros(size(X_t));
    ay = zeros(size(Y_t));
    az = zeros(size(Z_t));

    kappa_ob = zeros(N, 1);        % 初始化曲率观测
    tau_ob = zeros(N, 1);          % 初始化挠率观测

    for k = 1:N
        v = [XV_t(k), YV_t(k), ZV_t(k)];     % 速度向量
        a = [XA_t(k), YA_t(k), ZA_t(k)];     % 加速度向量
        aax = [0; diff(XA_t)];
        aay = [0; diff(YA_t)];
        aaz = [0; diff(ZA_t)];
        j = [aax(k), aay(k), aaz(k)];  % 加加速度向量
        
        v_cross_a = cross(v, a);                % 计算速度矢量和加速度矢量的叉乘
        norm_v = norm(v);                       % 速度矢量归一化
        norm_v_cross_a = norm(v_cross_a);       % 叉乘结果归一化
        
        % 曲率
        if norm_v > 1e-9
            kappa_ob(k) = norm_v_cross_a / (norm_v^3);% 曲率计算公式
        else
            kappa_ob(k) = 0;
        end
        
        % 挠率，这里采用的是知乎的计算公式变体 https://www.zhihu.com/question/25673951
        if norm_v_cross_a > 1e-9
            tau_ob(k) = dot(v_cross_a, j) / (norm_v_cross_a^2);% 扰率计算公式变体
        else
            tau_ob(k) = 0;
        end
    end

    kappa_ob = kappa_ob + normrnd(0, sigma_r, size(XV_t));  % 加噪声，作为观测
    tau_ob = tau_ob + normrnd(0, sigma_r, size(XV_t));

    
    rmse_obs = rmse_obs + sqrt(mean((x(obs_idx_beg:end)-X_t(obs_idx_beg:end)).^2 + ...
                         (y(obs_idx_beg:end)-Y_t(obs_idx_beg:end)).^2 + ...
                         (z(obs_idx_beg:end)-Z_t(obs_idx_beg:end)).^2));

    % 初始化变量
    Z = [x';vx';ax';y';vy';ay';z';vz';az';kappa_ob';tau_ob'];      % 观测数据
    X_pre = X0;
    Pk = P0;
    kappa = 0;
    tau = 0;
    I = eye(size(P0));
    tic;

    for k = 1:N
        % 在线计算卡尔曼增益 以及 预测&滤波的误差协方差矩阵
        % 先计算kappa 和 tau

        if k >= 3
            v = [X_cv_filted(2,k-1), X_cv_filted(5,k-1), X_cv_filted(8,k-1)]; % 速度向量
            a = [(X_cv_filted(2,k-1) - X_cv_filted(2,k-2)) / dt, 
                (X_cv_filted(5,k-1) - X_cv_filted(5,k-2)) / dt,
                (X_cv_filted(8,k-1) - X_cv_filted(8,k-2)) / dt];  % 加速度向量
            j = [(X_cv_filted(3,k-1) - X_cv_filted(3,k-2)) / dt, 
                (X_cv_filted(6,k-1) - X_cv_filted(6,k-2)) / dt,
                (X_cv_filted(9,k-1) - X_cv_filted(9,k-2)) / dt];  % 加加速度向量
        else
            v = [X_pre(2), X_pre(5), X_pre(8)];     % 速度向量
            a = [X_pre(3), X_pre(6), X_pre(9)];     % 加速度向量
            j = [0,0,0];
        end
        [kappa, tau] = kappa_tau_cal(v,a,j);        % 计算得到此时的曲率和扰率，由于是cv模型，加速度始终是0，所以tau始终是0

        Pk_estm = F_cv_augmented * Pk * F_cv_augmented' + Q_cv_augmented;                % 预测
        K_cv = Pk_estm * H_cv_augmented' / (H_cv_augmented * Pk_estm * H_cv_augmented' + R_cv);   % 卡尔曼增益
        Pk = (I - K_cv * H_cv_augmented) * Pk_estm;                           % 更新

        X_pre = KalmanFilter(Z(:,k), X_pre, F_cv_augmented, H_cv_augmented, K_cv);
        X_cv_filted(:,k) = X_pre + X_cv_filted(:,k);
    end
    cost_time(num) = toc;
    % disp("第"+num2str(num)+"次蒙特卡洛实验结束.");
end
disp("【蒙特卡洛实验结束!】 共进行"+num2str(Num)+"次实验。");
X_cv_filted = X_cv_filted ./ Num;
rmse_obs = rmse_obs ./ Num;
cost_time = mean(cost_time);
disp("平均实验耗时：" + num2str(cost_time*1000) + " ms");

x_filted = X_cv_filted(1,:);
xv_filted = X_cv_filted(2,:);
y_filted = X_cv_filted(4,:);
yv_filted = X_cv_filted(5,:);
z_filted = X_cv_filted(7,:);
zv_filted = X_cv_filted(8,:);

% 位置估计RMSE
rmse_xyz_cv = sqrt(mean((x_filted(obs_idx_beg:end)' - X_t(obs_idx_beg:end)).^2 + ...
                        (y_filted(obs_idx_beg:end)' - Y_t(obs_idx_beg:end)).^2 + ...
                        (z_filted(obs_idx_beg:end)' - Z_t(obs_idx_beg:end)).^2));
disp("=== 观测数据平均位置RMSE："+num2str(rmse_obs)+" m");
disp("=== 卡尔曼滤波后位置RMSE："+num2str(rmse_xyz_cv)+" m");

%% CV模型可视化滤波结果
vis = false;
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
    title("CV with kappa and tau");
    xlabel("x(m)"); ylabel("y(m)"); zlabel("z(m)");
    grid on; hold on;
    plot3(x_filted, y_filted, z_filted, 'LineWidth', 2);
    plot3(X_t, Y_t, Z_t, 'green', 'LineWidth', 2);
    legend("Observed Data", "Result", "Real Trajectory");
    exportgraphics(gcf, "策略二：CV+kappa+tau滤波结果.pdf", "ContentType", "vector")
    

    figure('name',"CV模型的速度卡尔曼滤波结果",'Position', [1300, 100, 1200, 800]);
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
    exportgraphics(gcf, "策略二：CV+kappa+tau模型的速度的滤波结果.pdf", "ContentType", "vector")

    
    figure('name',"CV模型的位置卡尔曼滤波结果",'Position', [1300, 100, 1200, 800]);
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
    exportgraphics(gcf, "策略二：CV+kappa+tau模型的位置的滤波结果.pdf", "ContentType", "vector")



    disp("可视化结果完成.");
end