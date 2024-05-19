%脚本，用于数据绘图
clear all;
close;
clc;
%数据为：
% 100  0.98146
% 200  0.95728
% 300  0.93369
% 400  0.91068
% 500  0.88824
% 600  0.86635

%数据绘图
x = [100, 200, 300, 400, 500, 600];
y = [0.98146, 0.95728, 0.93369, 0.91068, 0.88824, 0.86635];
plot(x, y, 'r', 'LineWidth', 2);
xlabel('时间步数');
ylabel('波场幅度');
title('介质1中波场幅度随时间步数变化');
grid on;
%在图上显示数据点
hold on;
plot(x, y, 'bo', 'MarkerSize', 10);
%显示数据标签
for i = 1:length(x)
    text(x(i), y(i), num2str(y(i)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end

%计算平均斜率
slope = diff(y) ./ diff(x);
mean_slope = mean(slope);
fprintf('平均斜率为：%f\n', mean_slope);

%第二组数据

% 800  0.55854
% 900  0.53187
% 1000  0.51454

%数据绘图
x = [ 800, 900, 1000];
y = [ 0.55854, 0.53187, 0.51454];
figure;
plot(x, y, 'r', 'LineWidth', 2);
xlabel('时间步数');
ylabel('波场幅度');
title('介质2中波场幅度随时间步数变化');
grid on;
%在图上显示数据点
hold on;
plot(x, y, 'bo', 'MarkerSize', 10);
%显示数据标签
for i = 1:length(x)
    text(x(i), y(i), num2str(y(i)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end

%计算平均斜率
slope = diff(y) ./ diff(x);
mean_slope = mean(slope);
fprintf('平均斜率为：%f\n', mean_slope);