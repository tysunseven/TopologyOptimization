function [mu_x_all, mu_y_all, path_distance] = generate_band_path_triangle()
    % 第一段 Γ→X (0,0)->(pi,0)
    x1 = linspace(0, pi, 20);
    y1 = zeros(1, 20);
    d1 = linspace(0, pi, 20);

    % 第二段 X→M (pi,0)->(pi,pi)
    x2 = pi * ones(1, 20);
    y2 = linspace(0, pi, 20);
    d2 = linspace(pi, 2*pi, 20);

    % 第三段 M→Γ (pi,pi)->(0,0) 直线
    x3 = linspace(pi, 0, 28);
    y3 = linspace(pi, 0, 28);
    seg3_length = hypot(pi, pi); % 直线距离为pi*sqrt(2)
    d3 = linspace(2*pi, 2*pi + seg3_length, 28);

    % 合并所有数据
    mu_x_all = [x1, x2, x3];
    mu_y_all = [y1, y2, y3];
    path_distance = [d1, d2, d3];

    % 归一化路径参数
    path_distance = (path_distance - min(path_distance)) / ...
                   (max(path_distance) - min(path_distance));
end