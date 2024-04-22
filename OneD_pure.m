import math.*

% 一维波动方程
% ∂^2 P / ∂ t^2 = c^2 * ∂^2 P / ∂ z^2

% 中心差分法离散化波动方程
% ∂^2 P / ∂ t^2 ≈ (P(t+Δt,z) - 2*P(t,z) + P(t-Δt,z)) / Δt^2
% ∂^2 P / ∂ z^2 ≈ (P(t,z+Δz) - 2*P(t,z) + P(t,z-Δz)) / Δz^2

% 代入一维波动方程，得到
% P(t+Δt,z) = 2*P(t,z) - P(t-Δt,z) + (c*Δt/Δz)^2 * (P(t,z+Δz) - 2*P(t,z) + P(t,z-Δz))

%CFL条件：c*Δt/Δz <= 1
% 设定参数
wave_speed = 1; % 波速
space_length = 10; % 空间长度
time_length = 10; % 时间长度
space_grid_num = 1000; % 空间网格数
time_grid_num = 1000; % 时间网格数
dz = space_length / space_grid_num; % 空间步长
dt = time_length / time_grid_num; % 时间步长
P = zeros(space_grid_num, time_grid_num); % 零初始化P域

% 对于平面正弦波： p = p0 * sin(ωt - kz)
% 色散公式： k^2 = ρ0 * K * ω^2  =>  k = ω / c

p0 = 1; % 振幅
omega = 2*pi; % 角频率
k = omega / wave_speed; % 波数

% O(h^2)精度的中心差分法，迭代计算
for m = 2:time_grid_num-1 % 时间迭代
    if m * dt < 1
        P(1, m) = p0 * sin(omega * m * dt); % 左边界（波源）
    end
    for n = 2:space_grid_num-1 % 空间迭代
        % 边界条件
        P(600, :) = 0; % z=6设定右边界
        P(n, m+1) = 2 * P(n, m) - P(n, m-1) + (wave_speed * dt / dz)^2 * (P(n+1, m) - 2 * P(n, m) + P(n-1, m));
        
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








