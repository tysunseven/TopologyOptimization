function figHandles = initOptimizationFigures(x, num_modes)
% 初始化所有图形组件并返回图形句柄结构体

%% 创建主图形窗口
fig = figure('Position', [0 170 1300 1200], 'Color', 'w');

%% ======================== 子图1：能带结构 ========================
ax_band = subplot(2,3,[1 2]);
colors = lines(num_modes);
hold(ax_band, 'on');
bandLines = gobjects(num_modes, 1);
for l = 1:num_modes
    bandLines(l) = plot(nan, nan, 'Color', colors(l,:), 'LineWidth', 1.5);
end

band_plots = gobjects(num_modes, 1);
for l = 1:num_modes
    band_plots(l) = plot(ax_band, NaN, NaN,...
                        'Color', colors(l,:),...
                        'LineWidth', 1.5);
end

% 带隙虚线初始化
gapLines = gobjects(2,1);
gapLines(1) = plot(nan, nan, '--', 'Color', [0.7 0.2 0.2], 'LineWidth', 1.5);
gapLines(2) = plot(nan, nan, '--', 'Color', [0.7 0.2 0.2], 'LineWidth', 1.5);

% 坐标轴设置
set(gca, 'XTick', [0, 0.292893, 0.585786, 1],...
         'XTickLabel', {'Γ','X','M','Γ'},...
         'XLim', [0 1],...
         'YLim', [0 7],...
         'FontSize', 12,...
         'FontWeight', 'bold');
ylabel('Frequency (Hz)');
title('Brillouin Zone Band Structure');
grid on;

%% 预创建带隙的Patch对象
band_gap_patches = gobjects(num_modes-1, 1); % 每个带隙对应一个Patch
for i = 1:num_modes-1
        % 初始创建不可见的Patch，覆盖整个X轴范围
    band_gap_patches(i) = patch(ax_band, 'XData', [0; 1; 1; 0], 'YData', [0; 0; 0; 0], ...
        'FaceColor', [1, 0.8, 0.8], 'EdgeColor', 'none', ...
        'FaceAlpha', 0.3, 'Visible', 'off');
    uistack(band_gap_patches(i), 'bottom'); % 确保在能带曲线下方
end

%% ======================== 子图2：密度图 ========================
subplot(2,3,3);
densityImage = imagesc(real(1-x));
colormap(gray);
caxis([0 1]);
axis equal tight;
axis off;
title('Density Plot');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');

%% ======================== 子图3：参数演化 ========================
ax3 = subplot(2,3,[4 5 6]);
hold(ax3, 'on');

% 定义图形样式
color_group1 = [0.2 0.6 0.9]; % 蓝色系
color_group2 = [0.9 0.3 0.2]; % 红色系
marker_size = 8;

% 创建绘图对象
h_max = plot(ax3, nan, nan, 'o-', 'Color', color_group1, 'MarkerSize', marker_size);
h_approxmax = plot(ax3, nan, nan, 's--', 'Color', color_group1, 'MarkerSize', marker_size);
h_min = plot(ax3, nan, nan, 'o-', 'Color', color_group2, 'MarkerSize', marker_size);
h_approxmin = plot(ax3, nan, nan, 's--', 'Color', color_group2, 'MarkerSize', marker_size);
h_approxgap = plot(ax3, nan, nan, '^-.', 'Color', [0.4 0.8 0.3], 'MarkerSize', marker_size);
h_realgap = plot(ax3, nan, nan, 'v:', 'Color', [0.4 0.8 0.3], 'MarkerSize', marker_size);

% 坐标轴设置
xlabel(ax3, '迭代次数', 'FontWeight','bold');
ylabel(ax3, '参数值', 'FontWeight','bold');
title(ax3, '优化参数演化过程', 'FontSize',12);
grid(ax3, 'on');
legend(ax3, {'Max(ω₁)','ApproxMax(ω₁)','Min(ω₂)','ApproxMin(ω₂)','估计带隙','实际带隙'},...
       'Location','bestoutside');

%% 封装图形句柄
figHandles = struct(...
    'fig', fig,...
    'ax_band', ax_band,...
    'bandLines', bandLines,...
    'gapLines', gapLines,...
    'band_gap_patches', band_gap_patches,...
    'band_plots', band_plots,...
    'densityImage', densityImage,...
    'ax3', ax3,...
    'h_max', h_max,...
    'h_approxmax', h_approxmax,...
    'h_min', h_min,...
    'h_approxmin', h_approxmin,...
    'h_approxgap', h_approxgap,...
    'h_realgap', h_realgap);
end