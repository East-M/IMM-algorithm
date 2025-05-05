function [xk] = KalmanFilter(zk, x_pre, Fk, Hk, K)
%KALMANFILTER 卡尔曼滤波器
%   对输入的数据zk进行卡尔曼滤波。
%   输入当前时刻的观测数据zk、上一时刻的估计值、状态转移矩阵Fk、量测矩阵Hk和离线计算好的卡尔曼增益K。
%   输出当前时刻的观测数据zk的滤波结果xk。
%
%示例：
%   [xk] = KalmanFilter(zk, x_pre, Fk, Hk, K)
    
    xk_estm = Fk * x_pre;                    % 预测
    xk = xk_estm + K * (zk - Hk * xk_estm);  % 更新
end

