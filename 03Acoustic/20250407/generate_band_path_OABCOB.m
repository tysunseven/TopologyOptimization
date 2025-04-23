function [mu_x_all, mu_y_all, path_distance] = generate_band_path_OABCOB()
    % 各段点数配置
    n_segment = 20;  % 前四段每段点数
    n_diag = 28;     % 第五段对角线点数
    
    % ================== 第一段 O→A (0,0)->(π,0) ==================
    x1 = linspace(0, pi, n_segment);
    y1 = zeros(1, n_segment);
    d1 = linspace(0, pi, n_segment); % 路径长度π
    
    % ================== 第二段 A→B (π,0)->(π,π) ==================
    x2 = pi * ones(1, n_segment);
    y2 = linspace(0, pi, n_segment);
    d2 = linspace(pi, 2*pi, n_segment); % 累计长度2π
    
    % ================== 第三段 B→C (π,π)->(0,π) ==================
    x3 = linspace(pi, 0, n_segment);
    y3 = pi * ones(1, n_segment);
    d3 = linspace(2*pi, 3*pi, n_segment); % 累计长度3π
    
    % ================== 第四段 C→O (0,π)->(0,0) ==================
    x4 = zeros(1, n_segment);
    y4 = linspace(pi, 0, n_segment);
    d4 = linspace(3*pi, 4*pi, n_segment); % 累计长度4π
    
    % ================== 第五段 O→B (0,0)->(π,π) ==================
    x5 = linspace(0, pi, n_diag);
    y5 = linspace(0, pi, n_diag);
    diag_length = pi*sqrt(2); % 对角线长度
    d5 = linspace(4*pi, 4*pi + diag_length, n_diag);
    
    % ================== 合并所有路径 ==================
    mu_x_all = [x1, x2, x3, x4, x5];
    mu_y_all = [y1, y2, y3, y4, y5];
    path_distance = [d1, d2, d3, d4, d5];
    
    % ================== 归一化路径参数 ==================
    path_distance = (path_distance - min(path_distance)) / ...
                   (max(path_distance) - min(path_distance));
end