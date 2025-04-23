function [figHandles, max_vals, min_vals, approxmax_vals, approxmin_vals, approxgap_vals, gap_vals] = ...
        updateFigures(figHandles, eigenvalues, x, mode, loop, ...
                                  path_distance, boundary_mu_x, num_modes, ...
                                  approxmax, approxmin, max_vals, min_vals, ...
                                  approxmax_vals, approxmin_vals, approxgap_vals, gap_vals)
% UPDATEOPTIMIZATIONFIGURES 更新优化过程的可视化图形
% 输入参数：
%   figHandles       - 图形句柄结构体
%   eigenvalues      - 特征值矩阵 (num_kpoints × num_modes)
%   x                - 设计变量矩阵
%   mode             - 目标模式编号
%   loop             - 当前迭代次数
%   path_distance    - 波矢路径坐标
%   boundary_mu_x    - 边界点x坐标
%   num_modes        - 总模式数
%   approxmax        - 近似最大值
%   approxmin        - 近似最小值
%   max_vals         - 历史最大值数组
%   min_vals         - 历史最小值数组
%   approxmax_vals   - 历史近似最大值数组
%   approxmin_vals   - 历史近似最小值数组
%   approxgap_vals   - 历史近似带隙数组
%   gap_vals         - 历史实际带隙数组
%
% 输出参数：
%   更新后的所有历史数据数组

%% 更新能带曲线
for l = 1:num_modes
    set(figHandles.bandLines(l),...
        'XData', path_distance,...
        'YData', eigenvalues(1:numel(boundary_mu_x), l)');
end

%% 更新带隙填充区域
for i = 1:num_modes-1
    % 提取当前模式对的特征值
    omega_i = eigenvalues(1:numel(boundary_mu_x), i);
    omega_i_plus_1 = eigenvalues(1:numel(boundary_mu_x), i+1);
    
    % 计算极值
    current_max = max(omega_i);
    current_min = min(omega_i_plus_1);
    
    % 更新图形对象
    if current_max < current_min
        set(figHandles.band_gap_patches(i),...
            'YData', [current_max, current_max, current_min, current_min],...
            'Visible', 'on');
    else
        set(figHandles.band_gap_patches(i), 'Visible', 'off');
    end
end

%% 更新密度分布图
set(figHandles.densityImage,...
    'CData', real(1 - x)); % 密度与设计变量反向显示

%% 记录历史数据
% 更新存储数组
max_vals(end+1) = max(eigenvalues(:, mode));
min_vals(end+1) = min(eigenvalues(:, mode+1));
approxmax_vals(end+1) = approxmax;
approxmin_vals(end+1) = approxmin;
approxgap_vals(end+1) = approxmax - approxmin;
gap_vals(end+1) = min(eigenvalues(:, mode+1)) - max(eigenvalues(:, mode));

%% 更新演化曲线
iterations = 1:loop;
set(figHandles.h_max,...
    'XData', iterations,...
    'YData', max_vals);
set(figHandles.h_approxmax,...
    'XData', iterations,...
    'YData', approxmax_vals);
set(figHandles.h_min,...
    'XData', iterations,...
    'YData', min_vals);
set(figHandles.h_approxmin,...
    'XData', iterations,...
    'YData', approxmin_vals);
set(figHandles.h_approxgap,...
    'XData', iterations,...
    'YData', approxgap_vals);
set(figHandles.h_realgap,...
    'XData', iterations,...
    'YData', gap_vals);

% 调整坐标轴范围
xlim(figHandles.ax3, [1, max(iterations)+0.5]);

%% 强制图形刷新
drawnow limitrate;
end


% for l = 1:num_modes
    %     set(figHandles.bandLines(l),...
    %         'XData', path_distance,...
    %         'YData', eigenvalues(1:numel(boundary_mu_x), l)');
    % end
    % 
    % for i = 1:num_modes-1
    %     % 提取当前模式i和i+1在所有k点中的频率
    %     omega_i = eigenvalues(1:numel(boundary_mu_x), i);
    %     omega_i_plus_1 = eigenvalues(1:numel(boundary_mu_x), i+1);
    % 
    %     % 计算全局最大值和最小值
    %     max_omega_i = max(omega_i);
    %     min_omega_i_p1 = min(omega_i_plus_1);
    % 
    %     % 判断是否存在带隙
    %     if max_omega_i < min_omega_i_p1
    %         % 设置Patch的Y坐标并显示
    %         set(figHandles.band_gap_patches(i), 'YData', [max_omega_i, max_omega_i, min_omega_i_p1, min_omega_i_p1], 'Visible', 'on');
    %     else
    %         set(figHandles.band_gap_patches(i), 'Visible', 'off');
    %     end
    % end
    % 
    % 
    % % 更新密度图
    % set(figHandles.densityImage, 'CData', real(1-x));
    % 
    % % 更新参数演化图
    % max_vals(end+1) = max(eigenvalues(:,mode));
    % approxmax_vals(end+1) = approxmax;
    % min_vals(end+1) = min(eigenvalues(:,mode+1));
    % approxmin_vals(end+1) = approxmin;
    % approxgap_vals(end+1) = approxmax - approxmin;
    % gap_vals(end+1) = min(eigenvalues(:,mode+1)) - max(eigenvalues(:,mode));
    % 
    % iterations = 1:loop;
    % set(figHandles.h_max, 'XData',iterations, 'YData',max_vals);
    % set(figHandles.h_approxmax, 'XData',iterations, 'YData',approxmax_vals);
    % set(figHandles.h_min, 'XData',iterations, 'YData',min_vals);
    % set(figHandles.h_approxmin, 'XData',iterations, 'YData',approxmin_vals);
    % set(figHandles.h_approxgap, 'XData',iterations, 'YData',approxgap_vals);
    % set(figHandles.h_realgap, 'XData',iterations, 'YData',gap_vals);
    % xlim(figHandles.ax3, [1 max(iterations)+0.5]);
    % 
    % drawnow limitrate;