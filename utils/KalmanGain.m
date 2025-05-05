function [ K ] = KalmanGain(P0, Fk, Hk, Qk, Rk)
%KALMANGAINOFFLINE 计算卡尔曼增益
%   根据输入的初始后验估计协方差矩阵P0、状态转移矩阵Fk、量测矩阵Hk、
%   过程噪声协方差矩阵Qk和量测噪声协方差矩阵Rk来离线计算卡尔曼增益K。
%   使用前需要调整该函数定义内输出卡尔曼增益的长度len。
%示例：
%   [K] = KalmanGain(P0, Fk, Hk, Qk, Rk)
    
    % 输出的卡尔曼增益长度
    len = 389;

    % 初始化变量
    tmp = P0 * Hk' / (Hk * P0 * Hk' + Rk);
    K = zeros([size(tmp), len]);       % 卡尔曼增益
    Pk = P0;
    I = eye(size(P0));

    % 离线计算卡尔曼增益
    for k = 1:len
        Pk_estm = Fk * Pk * Fk' + Qk;       % 预测
        K(:,:,k) = Pk_estm * Hk' / (Hk * Pk_estm * Hk' + Rk); % 卡尔曼增益
        Pk = (I - K(:,:,k) * Hk) * Pk_estm; % 更新
    end
end

