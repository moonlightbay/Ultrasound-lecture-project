%此代码在一维 P 波波动方程的基础上，实现了波动的反射/透射现象（无衰减）。
import math.*
clear all;
%界面反射系数
%R = p+ / p- = (Z2 - Z1) / (Z2 + Z1)
%T = 2 * Z2 / (Z2 + Z1)

%声学边界条件
%Z = ρ * c
%p0_1 = p0_2
%u1 = u2
%P1 = P2

%平面声波垂直入射到界面上，反射和透射的波动方程
%P+ = p0 * sin(ωt - k1z)      入射波
%P- = R * p0 * sin(ωt + k1z)   反射波
%Pt+ = T * p0 * sin(ωt - k2z)   透射波

% 设定基本参数
space_length = 10; % 空间长度
time_length = 10; % 时间长度
space_grid_num = 1000; % 空间网格数
time_grid_num = 2000; % 时间网格数
dz = space_length / space_grid_num; % 空间步长
dt = time_length / time_grid_num; % 时间步长
P = zeros(space_grid_num, time_grid_num); % 零初始化P域

% 设定介质参数
rho1 = 1; % 密度1
c1 = 1; % 波速1
rho2 = 2; % 密度2
c2 = 1.2; % 波速2
Z1 = rho1 * c1; % 阻抗1
Z2 = rho2 * c2; % 阻抗2
R = (Z2 - Z1) / (Z2 + Z1); % 反射系数
T = 2 * Z2 / (Z2 + Z1); % 透射系数

%反射系数与透射系数的关系
%R + T = 1

% 对于平面正弦波： p = p0 * sin(ωt - kz)
% 色散公式： k^2 = ρ0 * K * ω^2  =>  k = ω / c

p0 = 1; % 振幅
omega = 2 * pi; % 角频率
k1 = omega / c1; % 介质1波数
k2 = omega / c2; % 介质2波数


% O(h^2)精度的中心差分法，迭代计算
for m = 2:time_grid_num-1 % 时间迭代
    % 左边界处设置波源
    if m * dt < 1
    P(1, m) = p0 * sin(omega * m * dt);  % 简化波源模型为 p0 * sin(ωt)
    end
    for n = 2:space_grid_num-1 % 空间迭代
        if n < 600  % 介质1内
            P(n, m+1) = 2 * P(n, m) - P(n, m-1) + (c1 * dt / dz)^2 * (P(n+1, m) - 2 * P(n, m) + P(n-1, m));       
            pi = P(599,m+1);
            P(600, m+1) = R * pi * sin(omega*m*dt + k1*n*dz);
            P(601, m+1) = T * pi * sin(omega*m*dt - k2*n*dz);
        elseif n > 601  % 介质2内
            P(n, m+1) = 2 * P(n, m) - P(n, m-1) + (c2 * dt / dz)^2 * (P(n+1, m) - 2 * P(n, m) + P(n-1, m)); 
        end
    end
end



% 绘制波动图
figure;
for m = 1:time_grid_num
    plot(P(:, m));
    ylim([-5, 5]);
    line([600, 600], [-5, 5], 'Color', 'r', 'LineStyle', '--'); % 在x=600处画一条线
    drawnow;
    if strcmpi(get(gcf, 'currentkey'), 'q') % 按下q键退出
        break;
    end
end

% Close the figure and exit
close(gcf);