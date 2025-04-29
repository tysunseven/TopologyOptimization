function OutOfPlaneElasticBand(x)
nelx=size(x,1);
nely=size(x,2);
Lx=2;
h=Lx/nelx;
[space1,space2,time1,time2] = outofplaneconst();

num_modes = 10; num_random_points = 0;

kpoint = struct( 'mu_x', 0, 'mu_y', 0 );

[boundary_mu_x, boundary_mu_y, path_distance] = generate_band_path_OABCOB();

kpoints = repmat(kpoint, numel(boundary_mu_x)+num_random_points,1);
for i = 1:numel(boundary_mu_x)
    kpoints(i).mu_x = boundary_mu_x(i); kpoints(i).mu_y = boundary_mu_y(i);
end

eigenvalues = zeros(length(boundary_mu_x)+num_random_points, num_modes);
eigenvectors = zeros(length(boundary_mu_x)+num_random_points, (nelx+1)*(nely+1), num_modes);

[edofMat, Ke, Me, iIndex, jIndex] = init_fem(nelx,nely,h); [row, col, fixT] = init_trans(nelx,nely);
[fig, ax_band, ax_density, band_plots, density_plot, band_gap_patches] = init_plots(x, num_modes);

[K, M] = loop_fem(space1,space2,time1,time2,x,iIndex,jIndex,Ke,Me,nelx,nely);

parfor i = 1:length(kpoints)
    T = create_T(kpoints(i).mu_x, kpoints(i).mu_y, nelx, nely, row, col, fixT);
    K_tilde = T' * K * T; M_tilde = T' * M * T;
    [V, D] = eigs(K_tilde, M_tilde, num_modes, 'sm');
    eigenvalues(i,:) = sort(sqrt(abs(real(diag(D)))));
end

p = 48;                                 % 示例 p 值，根据实际需求调整
f = sum(eigenvalues(:, 1).^p).^(1/p);   % 手动计算 p-范数
disp(max(eigenvalues(:, 1)));
disp(f);


for l = 1:num_modes
    set(band_plots(l), 'XData', path_distance, 'YData', eigenvalues(1:numel(boundary_mu_x), l));
end

for i = 1:num_modes-1
    % 提取当前模式i和i+1在所有k点中的频率
    omega_i = eigenvalues(1:numel(boundary_mu_x), i);
    omega_i_plus_1 = eigenvalues(1:numel(boundary_mu_x), i+1);
    
    % 计算全局最大值和最小值
    max_omega_i = max(omega_i);
    min_omega_i_p1 = min(omega_i_plus_1);
    
    % 判断是否存在带隙
    if max_omega_i < min_omega_i_p1
        % 设置Patch的Y坐标并显示
        set(band_gap_patches(i), 'YData', [max_omega_i, max_omega_i, min_omega_i_p1, min_omega_i_p1], 'Visible', 'on');
    else
        set(band_gap_patches(i), 'Visible', 'off');
    end
end
    
set(density_plot, 'CData', real(1 - x)); drawnow limitrate;

