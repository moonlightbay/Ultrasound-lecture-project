%--------------------------------------------------------------------------
% 模拟二维声波在均匀介质中的传播
%--------------------------------------------------------------------------
% 参数设定
close all;
clear all;
clc;
Nx = 200;   % x方向的网格点数
Ny = 200;   % y方向的网格点数
c = 1500;   % 声速(m/s)
dt = 1e-6;  % 时间步长(s)
Nt = 2000;   % 时间步数
dx = 0.01;  % 空间步长(m)
dy = 0.01;  % 空间步长(m)

% 初始化
p = zeros(Nx, Ny);
p_prev = zeros(Nx, Ny);
p_next = zeros(Nx, Ny);

% 初始条件（波源）
x0 = 2;
% 阵列，间距为10
y0 = 95:1:105;
p(x0, y0) = 1;

% 时间演化
for t = 1:Nt
    for i = 2:Nx-1
        for j = 2:Ny-1
            % 使用二维波动方程的离散化公式           
            p_next(i, j) = 2*p(i, j) - p_prev(i, j) + c^2 * dt^2 / dx^2 * (p(i+1, j) + p(i-1, j) + p(i, j+1) + p(i, j-1) - 4*p(i, j));
        end
    end
    % 边界条件， 一阶吸收边界条件
    p_next(Nx, :) = p_next(Nx-1, :) - 1/c * (p_next(Nx-1, :) - p(Nx-1, :)) * dx/dt; %下边界
    p_next(:, Ny) = p_next(:, Ny-1) - 1/c * (p_next(:, Ny-1) - p(:, Ny-1)) * dy/dt; %右边界
    p_next(:, 1) = p_next(:, 2) - 1/c * (p_next(:, 2) - p(:, 2)) * dy/dt;

    % 更新波场
    p_prev = p;
    p = p_next;
    
    % 可视化
    imagesc(p);
    %加上坐标轴，标明xy方向
    xlabel('y');
    ylabel('x');
    colorbar;
    clim([-0.5,0.5]);
    pause(0.01);
    %按右上角退出绘图
    if ~ishandle(1)
        break;
    end
end

