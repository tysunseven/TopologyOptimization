function [max_mode1, min_mode2, max_mu_x_pi, max_mu_y_pi, min_mu_x_pi, min_mu_y_pi, gap]=postFigures(mode,num_modes,nelx,nely,row,col,fixT,K,M)

num_samples = 40; [mu_x_grid, mu_y_grid] = meshgrid(linspace(0, pi, num_samples), linspace(0, pi, num_samples));
mu_x_flat = mu_x_grid(:); mu_y_flat = mu_y_grid(:); eigenvalues_surface = zeros(num_samples^2, num_modes);

parfor i = 1:numel(mu_x_flat)
    T = create_T(mu_x_flat(i), mu_y_flat(i), nelx, nely, row, col, fixT);
    K_tilde = T' * K * T; M_tilde = T' * M * T;
    [~, D] = eigs(K_tilde, M_tilde, num_modes, 'sm');
    eigenvalues_surface(i,:) = sort(sqrt(abs(real(diag(D))))); 
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

[max_mode1, max_mode1_idx] = max(eigenvalues_surface(:,mode));
[min_mode2, min_mode2_idx] = min(eigenvalues_surface(:,mode+1));
max_mu_x_pi = mu_x_flat(max_mode1_idx)/pi; max_mu_y_pi = mu_y_flat(max_mode1_idx)/pi;
min_mu_x_pi = mu_x_flat(min_mode2_idx)/pi; min_mu_y_pi = mu_y_flat(min_mode2_idx)/pi;
gap = min_mode2 - max_mode1;

scatter3(max_mu_x_pi*pi, max_mu_y_pi*pi, max_mode1, 200, 'pentagram',... 
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','r',...
    'LineWidth',1.5);

scatter3(min_mu_x_pi*pi, min_mu_y_pi*pi, min_mode2, 200, 'o',...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'LineWidth',1.5);

fprintf('第%d级特征值全局最大值: %.4f，位置: (%.3fπ, %.3fπ)\n',...
        mode, max_mode1, max_mu_x_pi, max_mu_y_pi);
fprintf('第%d级特征值全局最小值: %.4f，位置: (%.3fπ, %.3fπ)\n',...
        mode+1, min_mode2, min_mu_x_pi, min_mu_y_pi);
fprintf('全局带隙宽度: %.4f\n\n', gap);