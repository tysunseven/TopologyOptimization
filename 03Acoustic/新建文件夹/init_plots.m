function [fig, ax_band, ax_density, band_plots, density_plot band_gap_patches] = init_plots(x, num_modes)
% 初始化拓扑优化可视化系统
% 输入参数:
%   x          - 初始设计变量矩阵 (nely × nelx)
%   num_modes  - 能带模式数量
% 输出参数:
%   fig          - 图形窗口句柄
%   ax_band      - 能带结构坐标轴句柄
%   ax_density   - 密度图坐标轴句柄
%   band_plots   - 能带曲线图形对象数组
%   density_plot - 密度图图像对象

%% 创建图形窗口
fig = figure('Position', [25 150 1400 600]);

%% ================= 能带结构子图初始化 =================
ax_band = subplot(1,2,1);
hold(ax_band, 'on');
colors = lines(num_modes);

% 预创建能带曲线对象
band_plots = gobjects(num_modes, 1);
for l = 1:num_modes
    band_plots(l) = plot(ax_band, NaN, NaN,...
                        'Color', colors(l,:),...
                        'LineWidth', 1.5);
end

% 配置坐标轴属性
set(ax_band, 'XTick', [0, 0.292893, 0.585786, 1],...
             'XTickLabel', {'Γ','X','M','Γ'},...
             'XLim', [0 1],...
             'YLim', [0 10],...
             'FontSize', 12,...
             'FontWeight', 'bold');
ylabel(ax_band, 'Frequency (Hz)');
title(ax_band, 'Brillouin Zone Band Structure');
grid(ax_band, 'on');

%% ================= 密度子图初始化 =================
ax_density = subplot(1,2,2);

% 创建初始密度图
density_plot = imagesc(ax_density, real(1 - x));
colormap(ax_density, gray);

% 配置坐标轴属性
set(ax_density, 'CLim', [0 1],...
                'DataAspectRatio', [1 1 1],...
                'Visible', 'off',...
                'FontSize', 12,...
                'FontWeight', 'bold');
title(ax_density, 'Density Plot');

%% 预创建带隙的Patch对象
band_gap_patches = gobjects(num_modes-1, 1); % 每个带隙对应一个Patch
for i = 1:num_modes-1
        % 初始创建不可见的Patch，覆盖整个X轴范围
    band_gap_patches(i) = patch(ax_band, 'XData', [0; 1; 1; 0], 'YData', [0; 0; 0; 0], ...
        'FaceColor', [1, 0.8, 0.8], 'EdgeColor', 'none', ...
        'FaceAlpha', 0.3, 'Visible', 'off');
    uistack(band_gap_patches(i), 'bottom'); % 确保在能带曲线下方
end






%% 强制首次绘制
drawnow;
end