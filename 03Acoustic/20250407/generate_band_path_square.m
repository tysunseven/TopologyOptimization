function [mu_x_all, mu_y_all, path_distance] = generate_band_path_square()
    % 第一段 Γ→X (0,0)->(pi,0)
    x1 = linspace(0, pi, 20);
    y1 = zeros(1, 20);
    d1 = linspace(0, pi, 20);

    % 第二段 X→M (pi,0)->(pi,pi)
    x2 = pi * ones(1, 20);
    y2 = linspace(0, pi, 20);
    d2 = linspace(pi, 2*pi, 20);

    % 第三段 M→Y (pi,pi)->(0,pi)
    x3 = linspace(pi, 0, 20);
    y3 = pi * ones(1, 20);
    d3 = linspace(2*pi, 3*pi, 20);

    % 第四段 Y→Γ (0,pi)->(0,0)
    x4 = zeros(1, 20);
    y4 = linspace(pi, 0, 20);
    d4 = linspace(3*pi, 4*pi, 20);

    % 合并所有数据
    mu_x_all = [x1, x2, x3, x4];
    mu_y_all = [y1, y2, y3, y4];
    path_distance = [d1, d2, d3, d4];

    % 归一化路径参数至 [0,1]
    path_distance = (path_distance - min(path_distance)) / ...
                   (max(path_distance) - min(path_distance));
end