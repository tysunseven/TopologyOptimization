function T = create_T(mu_x, mu_y, nelx, nely, row, col, fixT)
% 创建周期边界条件变换矩阵
% 输入参数:
%   mu_x  - X方向波矢
%   mu_y  - Y方向波矢
%   nelx  - X方向单元数
%   nely  - Y方向单元数
%   row   - 行索引集合
%   col   - 列索引集合
%   fixT  - 固定边界条件矩阵
% 输出参数:
%   T     - 周期边界变换矩阵

val = [repmat(exp(-1i*mu_y),1,nelx),exp(-1i*mu_x)*exp(-1i*mu_y),repmat(exp(-1i*mu_x),1,nely)];
T=sparse(row,col,val,(nely+1)*(nelx+1),(nely+1)*(nelx+1))+fixT;
T = T(:, ~ismember(1:size(T, 2), row));

end