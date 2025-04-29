function [x, nelx, nely, Lx, Ly, h] = init_x_from_png()
% 参数初始化函数
% 返回:
%   x     - 材料分布矩阵
%   nelx  - X方向单元数
%   nely  - Y方向单元数
%   Lx    - X方向尺寸
%   Ly    - Y方向尺寸
%   h     - 单元尺寸

%% 材料分布初始化
x_orig = imread("cellular1.jpg");      % 原始图像
x_orig = x_orig(:,:,1);                % 转换为灰度
x_orig = x_orig/255;                   % 归一化到[0,1]
x_orig = 1 - x_orig;                   % 反转颜色

% 自动裁剪边界
valid_rows = any(x_orig, 2);           % 检测有材料的行
valid_cols = any(x_orig, 1);           % 检测有材料的列
x_cropped = x_orig(valid_rows, valid_cols); 

% 转换为正方形域
[m, n] = size(x_cropped);
max_dim = max(m, n);                   % 取最大维度
pad_x = max_dim - m;                   % 竖直方向填充量
pad_y = max_dim - n;                   % 水平方向填充量

% 对称填充构建正方形矩阵
x_square = zeros(max_dim, max_dim, 'like', x_cropped);
x_square(1:m, 1:n) = x_cropped;        % 左上角填充原始数据
if pad_y > 0
    x_square(:, end-pad_y+1:end) = x_square(:, n:-1:1);  % 右侧对称填充
end
if pad_x > 0
    x_square(end-pad_x+1:end, :) = x_square(m:-1:1, :);  % 下侧对称填充
end

x = double(x_square);                 % 确保数据类型为double

%% 网格参数设置
Lx = 2;                               % X方向物理长度
Ly = Lx;                              % 保持与X方向一致
nelx = size(x, 1);                    % 单元数（基于正方形划分）
nely = nelx;                          % 保持与nelx一致
h = Lx / nelx;                        % 单元尺寸
end