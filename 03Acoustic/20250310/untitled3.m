%% 参数设置
A = 1;          % 振幅
k = 3*pi;       % 波数 
L = 1;          % 空间范围
omega = pi;   % 角频率
t_end = 20;      % 总时间
dt = 0.02;      % 时间步长
x = linspace(0, L, 200);  % 空间采样

%% 预定义波形函数
left_wave = @(t) A*cos(-k*x + omega*t);       % 左行波
right_wave = @(t) A*cos(k*x - 2*k*L + omega*t); % 右行波
combined_wave = @(t) left_wave(t) + right_wave(t); % 合成波

%% 图形界面初始化
figure('Color','white','Position',[100 100 800 900])
subplot(3,1,1)
h1 = plot(x, left_wave(0), 'b','LineWidth',1.5);
axis([0 L -A*1.2 A*1.2])
title('左行波: Acos(-kx + \omegat)')

subplot(3,1,2)
h2 = plot(x, right_wave(0), 'r','LineWidth',1.5);
axis([0 L -A*1.2 A*1.2])
title('右行波: Acos(kx - 2kL + \omegat)')

subplot(3,1,3)
h3 = plot(x, combined_wave(0), 'm','LineWidth',2);
axis([0 L -2*A*1.2 2*A*1.2])
title('合成驻波: Acos(kx-2kL+\omegat) + Acos(-kx+\omegat)')

%% 动态演示
for t = 0:dt:t_end
    % 计算各波形
    y1 = left_wave(t);
    y2 = right_wave(t);
    y3 = combined_wave(t);
    
    % 同步更新图形
    set(h1, 'YData', y1)
    set(h2, 'YData', y2)
    set(h3, 'YData', y3)
    
    % 添加统一时间显示
    subtitle(sprintf('波动传播演示  当前时间: %.2f s',t))
    drawnow
end