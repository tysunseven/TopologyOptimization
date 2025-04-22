function [space1,space2,time1,time2] = highcontrastconst()
% 计算材料的等效剪切模量和密度
% 输入：
%   x - 设计变量矩阵（0-1之间的值，表示材料分布）
% 输出：
%   mu - 等效剪切模量矩阵（与x同维度）
%   rho - 等效密度矩阵（与x同维度）

% 材料参数定义
rho1 = 1;       % 材料1密度
rho2 = 100;       % 材料2密度
E1 = 4;         % 材料1弹性模量
E2 = 400;        % 材料2弹性模量
nu = 0.34;      % 泊松比

% 计算剪切模量
mu1 = E1/(2*(1+nu));  % 材料1剪切模量
mu2 = E2/(2*(1+nu));  % 材料2剪切模量

space1 = mu1;
space2 = mu2;

time1 = rho1;
time2 = rho2;
end