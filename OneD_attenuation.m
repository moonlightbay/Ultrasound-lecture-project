% 此代码在一维P波波动方程（无衰减）实现波的反射/透射现象基础之上加入了声衰减效应。
% 按q可退出当前的图像。
% 衰减既不能只依赖空间步长，因为难以解决方向问题；也不能依赖于时间步长，因为难以知道什么时候衰减系数发生变化。
% P[n][m+1] = (2 * P[n][m] - P[n][m-1] * (1 - b * dt / 2) +(c * dt / dz) * (c * dt / dz) *(P[n+1][m] - 2 * P[n][m] + P[n-1][m])) /(1 + b * dt / 2);


close all;
clear all;
clc;

% 设定参数
space_length = 10; % 空间长度
time_length = 10; % 时间长度
space_grid_num = 1000; % 空间网格数
time_grid_num = 2000; % 时间网格数
dz = space_length / space_grid_num; % 空间步长
dt = time_length / time_grid_num; % 时间步长
P = zeros(space_grid_num, time_grid_num); % 零初始化P域
interface_position = 600; % 界面位置
p0 = 1; % 振幅
omega = 2*pi; % 角频率

% 设定声场参数
c1 = 1; % 介质1波速
c2 = 1.5; % 介质2波速
b = 0.05; % 空间衰减系数
alpha1 = 0.2;  % 介质1的空间衰减系数
alpha2 = 0.2;   % 介质2的空间衰减系数



% 迭代计算每个时间步
for m = 2:time_grid_num-1
    if m * dt < 0.5  % 左边界处设置波源，持续时间0.5秒
        P(1, m) = p0 * sin(omega * m * dt);  % 简化波源模型为 p0 * sin(ωt)
    end
    for n = 2:space_grid_num-1
        if n < interface_position  % 介质1内的波动方程，考虑衰减
            P(n, m+1) = (2 * P(n, m) - P(n, m-1) * (1 - b * dt / 2) +(c1 * dt / dz) * (c1 * dt / dz) *(P(n+1, m) - 2 * P(n, m) + P(n-1, m))) /(1 + b * dt / 2);
        else
             % 介质2内的波动方程，考虑衰减
            P(n, m+1) = (2 * P(n, m) - P(n, m-1) * (1 - b * dt / 2) +(c2 * dt / dz) * (c2 * dt / dz) *(P(n+1, m) - 2 * P(n, m) + P(n-1, m))) /(1 + b * dt / 2);
        end
    end   
    % 右边界,一阶吸收边界条件
    P(space_grid_num, m+1) = P(space_grid_num-1, m+1) - 1/c2 *(P(space_grid_num-1,m+1)-P(space_grid_num-1,m))*dz/dt;  
end

% 绘制波动图
figure;
for m = 1:3:time_grid_num % Increase the step size to plot at a higher frequency
    plot(P(:, m));
    ylim([-5, 5]);
    line([600, 600], [-5, 5], 'Color', 'r', 'LineStyle', '--'); % 在x=600处画一条线
    % 使用寻峰函数找到峰值并在图像上显示
    [pks, locs] = findpeaks(P(:, m),"MinPeakHeight",0.1);
    text(locs, pks+0.2, num2str(pks));
    
    drawnow;
    if strcmpi(get(gcf, 'currentkey'), 'q') % 按下q键退出
        break;
    end
end


% Close the figure and exit
close(gcf);




% 已废弃的方案，无法解决多方向衰减问题
% 不得不承认在多方向上的衰减是一个很复杂的问题，而我们用的只是一个简化的模型，可能无法仅通过增加一些项来解决这个问题。
% % 迭代计算每个时间步
% for m = 2:time_grid_num-1
%     if m * dt < 0.5  % 左边界处设置波源，持续时间0.5秒
%         P(1, m) = p0 * sin(omega * m * dt);  % 简化波源模型为 p0 * sin(ωt)
%     end
%     for n = 2:space_grid_num-1
%         if n < interface_position  % 介质1内的波动方程，考虑衰减
%             damping_factor = exp(-alpha1 * (n - 1) * dz);  % 计算衰减因子
%             P(n, m+1) = 2 * P(n, m) - P(n, m-1) + (c1 * dt / dz)^2 * damping_factor * (P(n+1, m) - 2 * P(n, m) + P(n-1, m));
%         else
%              % 介质2内的波动方程，考虑衰减
%             damping_factor = exp(-alpha2 * (n - interface_position) * dz);  % 计算衰减因子
%             P(n, m+1) = 2 * P(n, m) - P(n, m-1) + (c2 * dt / dz)^2 * damping_factor * (P(n+1, m) - 2 * P(n, m) + P(n-1, m));
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
%     % 使用寻峰函数找到峰值并在图像上显示
%     [pks, locs] = findpeaks(P(:, m),"MinPeakHeight",0.2);
%     text(locs+5, pks, num2str(pks));
    
%     drawnow;
%     if strcmpi(get(gcf, 'currentkey'), 'q') % 按下q键退出
%         break;
%     end
% end

% % Close the figure and exit
% close(gcf);


