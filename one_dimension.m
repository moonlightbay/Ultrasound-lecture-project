
import math.*

%中心差分法离散化波动方程
%∂^2 P / ∂ t^2 ≈ (P(t+Δt,z) - 2*P(t,z) + P(t-Δt,z)) / Δt^2
%∂^2 P / ∂ z^2 ≈ (P(t,z+Δz) - 2*P(t,z) + P(t,z-Δz)) / Δz^2

%代入波动方程，得到
%P(t+Δt,z) = 2*P(t,z) - P(t-Δt,z) + (c*Δt/Δz)^2 * (P(t,z+Δz) - 2*P(t,z) + P(t,z-Δz))

%设定参数
wave_speed = 1; %波速
space_length = 10; %空间长度
time_length = 5; %时间长度
space_grid_num = 1000; %空间网格数
time_grid_num = 1000; %时间网格数
dz = space_length/space_grid_num; %空间步长
dt = time_length/time_grid_num; %时间步长
P = zeros(space_grid_num, time_grid_num); %零初始化P域

%对于平面正弦波： p = p0 * sin(ωt - kz)
%色散公式： k^2 = ρ0 * K * ω^2  =>  k = ω / c
%初始化正弦波
p0 = 1; %振幅
omega = 200*pi; %角频率
k = omega / wave_speed; %波数
P(:,1) = p0 * sin(k*(1:space_grid_num)*dz); %初始化P域



%边界条件
P(1,:) = 0;
P(space_grid_num,:) = 0;


%O(h^2)精度的中心差分法，迭代计算
for m = 2:time_grid_num-1
    for n = 2:space_grid_num-1
        P(n,m+1) = 2*P(n,m) - P(n,m-1) + (wave_speed*dt/dz)^2 * (P(n+1,m) - 2*P(n,m) + P(n-1,m));
    end
end

%绘制波动图
figure;
for m = 1:time_grid_num
    plot(P(:,m));
    drawnow;
end








