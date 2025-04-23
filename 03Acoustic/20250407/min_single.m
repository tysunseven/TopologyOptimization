%% 将目标函数和灵敏度计算提取到单独文件（compute_objective.m）
function [f0val, df0dx, fval, dfdx] = min_single(numb, kpoints, eigenvalues, eigenvectors, mu2, mu1, rho2, rho1, edofMata, Ke, Me, M, x, n)
% 计算目标函数和灵敏度
% 输入参数说明：
%   boundary_points - 结构体数组，包含特征值和特征向量
%   mu2, mu1, rho2, rho1 - 材料参数 
%   edofMata, Ke, Me - 有限元模型参数
%   M - 全局质量矩阵
%   x - 设计变量矩阵
%   n - 设计变量总数

% 提取第40个边界点的基态信息
disp(kpoints(numb));
f = eigenvalues(numb,1);
phi = eigenvectors(numb,:,1)';

% 计算灵敏度
ceK = (mu2-mu1)*reshape(sum((phi(edofMata)*Ke).*phi(edofMata),2),size(x));
ceM = (rho2-rho1)*reshape(sum((phi(edofMata)*Me).*phi(edofMata),2),size(x));
dc = real((ceK-f*ceM)/(phi'*M*phi));

% 目标函数参数
f0val = f;                                    % 目标函数值
df0dx = dc(:);                                % 目标函数梯度
dv = ones(size(x)); 
fval  = -(sum(x(:))/(0.4*n) - 1);            % 体积约束函数值
dfdx  = -dv(:)'/ (0.4*n);                    % 约束梯度
end