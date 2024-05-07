% 此代码基于一维P波波动方程（无衰减）实现波的反射/透射现象。

% 一维波动方程
% ∂^2 P / ∂ t^2 = c^2 * ∂^2 P / ∂ z^2

% 使用中心差分法对波动方程进行离散化
% ∂^2 P / ∂ t^2 ≈ (P(t+Δt,z) - 2*P(t,z) + P(t-Δt,z)) / Δt^2
% ∂^2 P / ∂ z^2 ≈ (P(t,z+Δz) - 2*P(t,z) + P(t,z-Δz)) / Δz^2

% 将上述结果代入一维波动方程，得到
% P(t+Δt,z) = 2*P(t,z) - P(t-Δt,z) + (c*Δt/Δz)^2 * (P(t,z+Δz) - 2*P(t,z) + P(t,z-Δz))

% CFL条件：c*Δt/Δz <= 1

% 界面反射系数
% R = p+ / p- = (Z2 - Z1) / (Z2 + Z1)
% T = 2 * Z2 / (Z2 + Z1)

% 声阻抗相关公式
% Z = ρ * c = sqrt(ρ/K),    K为弹性模量，ρ为密度
% p0_1 = p0_2
% u1 = u2
% P1 = P2

% 垂直入射平面声波在界面上的反射和透射波动方程
% P+ = p0 * sin(ωt - k1z)      入射波
% P- = R * p0 * sin(ωt + k1z)   反射波
% Pt+ = T * p0 * sin(ωt - k2z)   透射波

% 对于平面正弦波：p = p0 * sin(ωt - kz)
% 色散公式：k^2 = ρ0 * K * ω^2  =>  k = ω / c

% 波速和密度的关系  c = sqrt(K / ρ),K为弹性模量，ρ为密度
% 这意味着波速由介质本身的性质决定

%一阶边界条件的离散化

%-----------------------------------------------------------------------------------

import math.*
clear all;
% 设定声场参数
space_length = 10; % 空间长度
time_length = 10; % 时间长度
space_grid_num = 1000; % 空间网格数
time_grid_num = 1000; % 时间网格数
dz = space_length / space_grid_num; % 空间步长
dt = time_length / time_grid_num; % 时间步长
P = zeros(space_grid_num, time_grid_num); % 零初始化P域
interface_position = 600; % 界面位置
c1 = 1; % 介质1波速
c2 = 0.8; % 介质2波速
p0 = 1; % 振幅
omega = 2 * pi; % 角频率

% k1 = omega / c1; % 介质1波数
% k2 = omega / c2; % 介质2波数
rho1 = 1; % 密度1
rho2 = 2; % 密度2
Z1 = rho1 * c1; % 阻抗1
Z2 = rho2 * c2; % 阻抗2
R = (Z2 - Z1) / (Z2 + Z1); % 反射系数
T = 2 * Z2 / (Z2 + Z1); % 透射系数

%--------------------------只改变了波速------------------------------------------
% %迭代计算每个时间步
% for m = 2:time_grid_num-1
%     if m * dt < 0.5  % 左边界处设置波源，持续时间0.5秒
%         P(1, m) = p0 * sin(omega * m * dt);  % 简化波源模型为 p0 * sin(ωt)
%     end
%     for n = 2:space_grid_num-1  %599
%         if n < interface_position  % 介质1内的波动方程
%             P(n, m+1) = 2 * P(n, m) - P(n, m-1) + (c1 * dt / dz)^2 * (P(n+1, m) - 2 * P(n, m) + P(n-1, m));
%         else
%              % 介质2内的波动方程
%             P(n, m+1) = 2 * P(n, m) - P(n, m-1) + (c2 * dt / dz)^2 * (P(n+1, m) - 2 * P(n, m) + P(n-1, m));
%         end
%     end   
%     % 右边界,一阶吸收边界条件
%     P(space_grid_num, m+1) = P(space_grid_num-1, m+1) - 1/c2 *(P(space_grid_num-1,m+1)-P(space_grid_num-1,m))*dz/dt;  
% end
% % 绘制波动图
% figure;
% for m = 1:3:time_grid_num % Increase the step size to plot at a higher frequency
%     plot(P(:, m));
%     ylim([-5, 5]);
%     line([600, 600], [-5, 5], 'Color', 'r', 'LineStyle', '--'); % 在x=600处画一条线
%     drawnow;
%     if strcmpi(get(gcf, 'currentkey'), 'q') % 按下q键退出
%         break;
%     end
% end
%-----------------------------------------------------------------------------------
%分别计算入射波、反射波和透射波
Pi = zeros(space_grid_num, time_grid_num); % 零初始化Pi域
Pt = zeros(space_grid_num, time_grid_num); % 零初始化Pt域
Pr = zeros(space_grid_num, time_grid_num); % 零初始化Pr域

%计算出到达界面的时间步
time_step = interface_position * dz / c1 / dt;

%计算入射波
for m = 2:time_step
    if m * dt < 0.5
        Pi(1, m) = p0 * sin(omega * m * dt); % 左边界（波源）
    end
    for n = 2:space_grid_num-1
        Pi(n, m+1) = 2 * Pi(n, m) - Pi(n, m-1) + (c1 * dt / dz)^2 * (Pi(n+1, m) - 2 * Pi(n, m) + Pi(n-1, m));
    end
    % 右边界,一阶吸收边界条件
    Pi(space_grid_num, m+1) = Pi(space_grid_num-1, m+1) - 1/c1 *(Pi(space_grid_num-1,m+1)-Pi(space_grid_num-1,m))*dz/dt;  
end
%绘制入射波
figure;
for m = 1:3:time_step % Increase the step size to plot at a higher frequency
    plot(Pi(:, m));
    ylim([-5, 5]);
    line([600, 600], [-5, 5], 'Color', 'r', 'LineStyle', '--'); % 在x=600处画一条线
    drawnow;
    if strcmpi(get(gcf, 'currentkey'), 'q') % 按下q键退出
        break;
    end
end



%计算反射波
for m = time_step:time_grid_num-1
    if m * dt < 0.5
        Pr(interface_position, m) = Pi(interface_position, m); % 反射波的初始值
    end
    for n = interface_position:space_grid_num-1
        Pr(n, m+1) = 2 * Pr(n, m) - Pr(n, m-1) + (c2 * dt / dz)^2 * (Pr(n+1, m) - 2 * Pr(n, m) + Pr(n-1, m));
    end
    % 右边界,一阶吸收边界条件
    Pr(space_grid_num, m+1) = Pr(space_grid_num-1, m+1) - 1/c2 *(Pr(space_grid_num-1,m+1)-Pr(space_grid_num-1,m))*dz/dt;  
end






% Close the figure and exit
close(gcf);


