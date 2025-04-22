function [x, nelx, nely, Lx, Ly, h] = init_x_from_zero()
% 参数初始化函数
% 返回:
%   x     - 材料分布矩阵
%   nelx  - X方向单元数
%   nely  - Y方向单元数
%   Lx    - X方向尺寸
%   Ly    - Y方向尺寸
%   h     - 单元尺寸

Lx = 2;                               % X方向物理长度
Ly = Lx;                              % 保持与X方向一致
nelx = 20;                    % 单元数（基于正方形划分）
nely = 20;                          % 保持与nelx一致
h = Lx / nelx;                        % 单元尺寸
x = zeros(nely, nelx);
end