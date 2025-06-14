%% ====================================================================== %
% Sun Yat-sen University 电子与通信工程学院
% 统计信号处理课程1班 - 第3小组
% @author: 成先锋 莫晓东 陈立邦 成泽宇
% @email:  moxd3@mail2.sysu.edu.cn
% @GitHub: https://github.com/East-M
% date:2025年5月5日
%
% #code: 大机动飞机目标的跟踪滤波技术性能仿真
% 解决问题二的利用IMM模型对数据进行滤波，并做可视化处理。
%% ====================================================================== %
clear;
clc;
close all;

addpath utils\ data\
%% NOTICE
disp("【大机动飞机目标的跟踪滤波技术性能仿真】");
disp("解决问题二的利用IMM模型对数据进行滤波，并做可视化处理。");
disp("选用的子模型为：CV模型  CA模型  Singer模型");
disp(" ");
%% 路径定义
K_cv_Path = "./data/K_cv_imm.mat";
K_ca_Path = "./data/K_ca_imm.mat";
K_singer_Path = "./data/K_singer_imm.mat";
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
X0 = [X_t(1); (X_t(2)-X_t(1))/T; (X_t(3)-2*X_t(2)+X_t(1))/T^2; % 初始位置、速度估计和加速度估计
          Y_t(1); (Y_t(2)-Y_t(1))/T; (Y_t(3)-2*Y_t(2)+Y_t(1))/T^2;
          Z_t(1); (Z_t(2)-Z_t(1))/T; (Z_t(3)-2*Z_t(2)+Z_t(1))/T^2];
P0 = diag([sigma_r^2, 50^2, 25^2, sigma_r^2, 50^2, 25^2, sigma_r^2, 50^2, 25^2]);

%% @CV模型参数
% 状态转移矩阵
F_cv = blkdiag([1,T,0;0,1,0;0,0,0],[1,T,0;0,1,0;0,0,0],[1,T,0;0,1,0;0,0,0]);

% 系统扰动矩阵
G_cv = blkdiag([0.5*T^2;T;0],[0.5*T^2;T;0],[0.5*T^2;T;0]);

% 量测矩阵
H_cv = blkdiag([1,0,0],[1,0,0],[1,0,0]);

% 过程噪声协方差矩阵
Q_cv = sigma^2 * (G_cv * G_cv');

% 量测噪声协方差矩阵
R_cv = sigma_r^2 * eye(zNum);

%% @CA模型参数
% 状态转移矩阵
F_ca = blkdiag([1,T,T^2/2;0,1,T;0,0,1],[1,T,T^2/2;0,1,T;0,0,1],[1,T,T^2/2;0,1,T;0,0,1]);

% 系统扰动矩阵
G_ca = blkdiag([T^3/6;T^2/2;T],[T^3/6;T^2/2;T],[T^3/6;T^2/2;T]);

% 量测矩阵
H_ca = blkdiag([1,0,0],[1,0,0],[1,0,0]);

% 过程噪声协方差矩阵
Q_ca = sigma^2 * (G_ca * G_ca');

% 量测噪声协方差矩阵
R_ca = sigma_r^2 * eye(zNum);

%% @Singer模型参数
sigma_m = 70;       % Singer模型的加速度方差
alpha = 1/5;        % 机动频率，经验值。

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

%% IMM模型的参数定义
% 模型个数
r = 3;

% 转移概率矩阵
P = [0.8,0.1,0.1;
     0.1,0.8,0.1;
     0.1,0.1,0.8];

% 初始模型概率
mu = [1/3;1/3;1/3];

% 子模型参数
% 状态转移矩阵
F(:,:,1) = F_cv;
F(:,:,2) = F_ca;
F(:,:,3) = F_singer;

% 过程噪声协方差矩阵
Q(:,:,1) = Q_cv;
Q(:,:,2) = Q_ca;
Q(:,:,3) = Q_singer;

% 量测矩阵
H(:,:,1) = H_cv;
H(:,:,2) = H_ca;
H(:,:,3) = H_singer;

% 量测噪声协方差矩阵
R(:,:,1) = R_cv;
R(:,:,2) = R_ca;
R(:,:,3) = R_singer;

%% IMM算法
X_IMM_filtered = zeros(xNum,N);
mu_record = zeros(r, N);
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

    % 初始值
    X_pre(:,1) = X0;
    X_pre(:,2) = X0;
    X_pre(:,3) = X0;
    
    P_pre(:,:,1) = P0;
    P_pre(:,:,2) = P0;
    P_pre(:,:,3) = P0;

    tic;
    for k = 1:N
        X_0 = zeros(xNum, r);
        P_0 = zeros(size(P0));
        v = zeros(zNum,r);
        S = zeros([zNum,zNum,r]);
        likelihood = zeros(r,1);

        % 输入交互
        for j = 1:r
            c = P(:,j)' * mu;
            X_0(:,j) = 1/c*(P(1,j) * mu(1) * X_pre(:,1) + ...
                            P(2,j) * mu(2) * X_pre(:,2) + ...
                            P(3,j) * mu(3) * X_pre(:,3));
            P_0(:,:,j) = 1/c*(P(1,j)*mu(1)*(P_pre(:,:,1)+(X_pre(:,1)-X_0(:,j))*(X_pre(:,1)-X_0(:,j))')+...
                P(2,j)*mu(2)*(P_pre(:,:,2)+(X_pre(:,2)-X_0(:,j))*(X_pre(:,2)-X_0(:,j))')+...
                P(3,j)*mu(3)*(P_pre(:,:,3)+(X_pre(:,3)-X_0(:,j))*(X_pre(:,3)-X_0(:,j))'));
        end

        % 模型条件滤波（对每个模型j独立进行卡尔曼滤波）
        for j = 1:r
            % 预测
            X_pre(:,j) = F(:,:,j) * X_0(:,j);
            P_pre(:,:,j) = F(:,:,j) * P_0(:,:,j) * F(:,:,j)' + Q(:,:,j);

            % 更新
            v(:,j) = Z(:,k)-H(:,:,j)*X_pre(:,j);
            S(:,:,j) = H(:,:,j)*P_pre(:,:,j)*H(:,:,j)'+R(:,:,j);

            K = P_pre(:,:,j) * H(:,:,j)' / S(:,:,j);
            X_pre(:,j) = X_pre(:,j) + K * v(:,j);
            P_pre(:,:,j) = (eye(xNum)-K*H(:,:,j)) * P_pre(:,:,j);
        end

        % 模型概率更新
        mu_tmp = zeros(size(mu));
        for j = 1:r
            likelihood(j) = 1 / sqrt(det(2*pi*S(:,:,j))) * exp(-1/2*(v(:,j)'/S(:,:,j))*v(:,j));
        end
        for j = 1:r
            mu_tmp(j) = likelihood(j)*P(:,j)'*mu / ...
                    (likelihood' * P' * mu);
        end
        mu = mu_tmp / sum(mu_tmp);
        
        % 输出融合结果
        X_IMM_filtered(:,k) = mu(1)*X_pre(:,1)+mu(2)*X_pre(:,2)+mu(3)*X_pre(:,3) + X_IMM_filtered(:,k);
        mu_record(:,k) = mu;
    end
    cost_time(num) = toc;
end
disp("【蒙特卡洛实验结束!】 共进行"+num2str(Num)+"次实验。");
X_IMM_filtered = X_IMM_filtered ./ Num;
rmse_obs = rmse_obs ./ Num;
cost_time = mean(cost_time);
disp("平均实验耗时：" + num2str(cost_time*1000) + " ms");

x_filted = X_IMM_filtered(1,:);
xv_filted = X_IMM_filtered(2,:);
y_filted = X_IMM_filtered(4,:);
yv_filted = X_IMM_filtered(5,:);
z_filted = X_IMM_filtered(7,:);
zv_filted = X_IMM_filtered(8,:);

% 位置估计RMSE
rmse_xyz_imm = sqrt(mean((x_filted(obs_idx_beg:end)' - X_t(obs_idx_beg:end)).^2 + ...
                        (y_filted(obs_idx_beg:end)' - Y_t(obs_idx_beg:end)).^2 + ...
                        (z_filted(obs_idx_beg:end)' - Z_t(obs_idx_beg:end)).^2));
disp("=== 观测数据平均位置RMSE："+num2str(rmse_obs)+" m");
disp("=== 卡尔曼滤波后位置RMSE："+num2str(rmse_xyz_imm)+" m");

%% IMM模型可视化滤波结果
figure('name',"实验数据",'Position', [100, 100, 1200, 800]);
plot3(X_t, Y_t, Z_t, 'LineWidth', 2); title("飞机真实轨迹");
xlabel("x方向 (m)"); ylabel("y方向 ((m)"); zlabel("z方向 (m)");
grid on; hold on;
plot3(x, y, z, '.', 'MarkerSize', 8); title("飞行轨迹观测数据");
% plot3(x_filted, y_filted, z_filted, 'LineWidth', 2);
legend("Real Trajectory", "Observed Data");
%%
figure('name',"滤波结果3维图",'Position', [100, 100, 1200, 800]);
plot3(x, y, z, '.', 'MarkerSize', 8);
title("the Result of IMM");
xlabel("x (m)"); ylabel("y ((m)"); zlabel("z (m)");
grid on; hold on;
plot3(x_filted, y_filted, z_filted, 'LineWidth', 2);
plot3(X_t, Y_t, Z_t, 'green', 'LineWidth', 2);
legend("Observed Data", "Result", "Real Trajectory");
exportgraphics(gcf, "IMM滤波结果.pdf", "ContentType", "vector")
%% 
figure('name',"IMM模型的速度卡尔曼滤波结果",'Position', [1300, 100, 1200, 800]);
subplot(3,1,1);
plot(XV_t); xlabel("k"); ylabel("m/s"); title("xv分量");
hold on; grid on;
plot(xv_filted);
legend("True Speed","Result");

subplot(3,1,2);
plot(YV_t); xlabel("k"); ylabel("m/s"); title("yv分量");
hold on; grid on;
plot(yv_filted);
legend("True Speed","Result");

subplot(3,1,3);
plot(ZV_t); xlabel("k"); ylabel("m/s"); title("zv分量");
hold on; grid on;
plot(zv_filted);
legend("True Speed","Result");
exportgraphics(gcf, "IMM速度的滤波结果.pdf", "ContentType", "vector")
%%
figure('name',"IMM模型的位置卡尔曼滤波结果",'Position', [1300, 100, 1200, 800]);
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
exportgraphics(gcf, "IMM位置的滤波结果.pdf", "ContentType", "vector")
%%
figure('name',"概率转化图", 'Position', [1300, 100, 1200, 800]);
plot(mu_record(1,:), 'LineWidth', 2); hold on; grid on;
plot(mu_record(2,:), 'LineWidth', 2);
plot(mu_record(3,:), 'LineWidth', 2);
xlabel("time"); ylabel("probability");
legend("CV", "CA", "Singer");
exportgraphics(gcf, "IMM滤波概率转化图.pdf", "ContentType", "vector")

disp("可视化结果完成.");