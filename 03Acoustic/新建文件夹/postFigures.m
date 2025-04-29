function [bdmax, bdmin, bdmax_x, bdmax_y, bdmin_x, bdmin_y, bdgap, gbmax, gbmin, gbmax_x, gbmax_y, gbmin_x, gbmin_y, gbgap]=postFigures(mode,num_modes,nelx,nely,row,col,fixT,K,M)

num_samples = 40; [mu_x_grid, mu_y_grid] = meshgrid(linspace(0, pi, num_samples), linspace(0, pi, num_samples));
mu_x_flat = mu_x_grid(:); mu_y_flat = mu_y_grid(:); eigenvalues_surface = zeros(num_samples^2, num_modes);

parfor i = 1:numel(mu_x_flat)
    T = create_T(mu_x_flat(i), mu_y_flat(i), nelx, nely, row, col, fixT);
    K_tilde = T' * K * T; M_tilde = T' * M * T;
    [~, D] = eigs(K_tilde, M_tilde, num_modes, 'sm');
    eigenvalues_surface(i,:) = sort(sqrt(abs(real(diag(D))))); 
end

%% 第二部分：从网格提取边界点（增强版）
epsilon = 1e-5; % 边界判断容差

% 创建扩展边界判断条件（新增x=y对角线）
is_boundary = (abs(mu_x_flat) < epsilon) | ...       % x=0
             (abs(mu_x_flat - pi) < epsilon) | ...  % x=π
             (abs(mu_y_flat) < epsilon) | ...       % y=0
             (abs(mu_y_flat - pi) < epsilon) | ...  % y=π
             (abs(mu_x_flat - mu_y_flat) < epsilon); % x=y对角线

% 提取边界数据
boundary_mu_x = mu_x_flat(is_boundary);
boundary_mu_y = mu_y_flat(is_boundary);
boundary_eigen = eigenvalues_surface(is_boundary, :);

%% 第三部分：边界特征分析
if ~isempty(boundary_eigen)
    % 当前模式的最大值
    [bdmax, max_idx] = max(boundary_eigen(:, mode));
    bdmax_x = boundary_mu_x(max_idx)/pi;
    bdmax_y = boundary_mu_y(max_idx)/pi;
    
    % 下一模式的最小值
    [bdmin, min_idx] = min(boundary_eigen(:, mode+1));
    bdmin_x = boundary_mu_x(min_idx)/pi;
    bdmin_y = boundary_mu_y(min_idx)/pi;
    bdgap = bdmin - bdmax;
else
    warning('未找到边界点，请检查网格采样设置');
    [bdmax, bdmin, bdgap] = deal(NaN);
    [bdmax_x, bdmax_y, bdmin_x, bdmin_y] = deal(NaN);
end

eigenvalues_grid = cell(num_modes, 1);
for l = 1:num_modes
    eigenvalues_grid{l} = reshape(eigenvalues_surface(:, l), num_samples, num_samples);
end

figure('Position', [1400 600 900 700]); colors = lines(num_modes); hold on;
for l = 1:num_modes
    surf(mu_x_grid, mu_y_grid, eigenvalues_grid{l}, 'FaceColor', colors(l,:), 'EdgeColor', 'none', 'FaceAlpha', 0.85);
end

colormap(lines(num_modes)); clim([1 num_modes]); lighting gouraud; camlight(30,30); view(-37.5, 30); axis tight; grid on;
xlabel('\mu_x (rad)', 'FontWeight', 'bold'); ylabel('\mu_y (rad)', 'FontWeight', 'bold'); zlabel('Frequency (Hz)', 'FontWeight', 'bold');

[gbmax, max_mode1_idx] = max(eigenvalues_surface(:,mode));
[gbmin, min_mode2_idx] = min(eigenvalues_surface(:,mode+1));
gbmax_x = mu_x_flat(max_mode1_idx)/pi; gbmax_y = mu_y_flat(max_mode1_idx)/pi;
gbmin_x = mu_x_flat(min_mode2_idx)/pi; gbmin_y = mu_y_flat(min_mode2_idx)/pi;
gbgap = gbmin - gbmax;

scatter3(gbmax_x*pi, gbmax_y*pi, gbmax, 200, 'pentagram',... 
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','r',...
    'LineWidth',1.5);

scatter3(gbmin_x*pi, gbmin_y*pi, gbmin, 200, 'o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'LineWidth',1.5);

fprintf('第%d级特征值全局最大值: %.4f，位置: (%.3fπ, %.3fπ)\n',...
        mode, gbmax, gbmax_x, gbmax_y);
fprintf('第%d级特征值全局最小值: %.4f，位置: (%.3fπ, %.3fπ)\n',...
        mode+1, gbmin, gbmin_x, gbmin_y);
fprintf('全局带隙宽度: %.4f\n\n', gbgap);