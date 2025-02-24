% 参数设置
nelx = 100;  % x方向网格数
nely = 50;  % y方向网格数
omega = 0.01;  % 频率
rho1 = 1;  % 空气的密度
kappa1 = 1;  % 空气的体积模量
rho2 = 1000;  % 材料的密度
kappa2 = 1000;  % 材料的体积模量
xi = zeros(nely, nelx);

% 找到中间列的索引
mid_col = nelx / 2;

% 找到中间两个元素的索引
mid_row1 = nely / 2;
mid_row2 = mid_row1 + 1;

% 在中间列的前 (nely/2)-1 个元素和后 (nely/2)-1 个元素赋值为 1
xi(1:mid_row1-1, mid_col) = 1;
xi(mid_row2+1:end, mid_col) = 1;

% 其他列保持为 0（默认值）
colormap(gray); imagesc(1-xi); caxis([0 1]); axis equal; axis off; drawnow;

% 材料属性
rho = rho2 * xi + rho1;
kappa = kappa2 * xi + kappa1;

%% PREPARE FINITE ELEMENT ANALYSIS
KE = 1/6*[4 -1 -2 -1;
     -1 4 -1 -2;
     -2 -1 4 -1;
     -1 -2 -1 4];

ME = 1/36*[4 2 1 2;
    2 4 2 1;
    1 2 4 2;
    2 1 2 4];

nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx); % 给(nely+1)*(nelx+1)个节点编号
edofVec = reshape(nodenrs(1:end-1,1:end-1)+nely+1,nelx*nely,1); % 记录每个单元的右上节点的编号
edofMat = repmat(edofVec,1,4)+repmat([0 -nely-1 -nely 1],nelx*nely,1); % 记录每个单元的四个节点的编号

row_nely = edofMat(nely, :);

% 显示结果
disp('第 nely 行的四个元素是:');
disp(row_nely);
iIndex = reshape(kron(edofMat,ones(4,1))',16*nelx*nely,1);
jIndex = reshape(kron(edofMat,ones(1,4))',16*nelx*nely,1);

sK = reshape(KE(:)*rho(:)',16*nelx*nely,1);
sM = reshape(ME(:)*kappa(:)',16*nelx*nely,1);

K = sparse(iIndex,jIndex,sK); K = (K+K')/2;
M = sparse(iIndex,jIndex,sM); M = (M+M')/2;

F = zeros((nely+1)*(nelx+1), 1);
F(1:nely+1) = 1;

disp('维度检查:');
disp(['K 的维度: ', num2str(size(K))]);
disp(['M 的维度: ', num2str(size(M))]);
disp(['F 的维度: ', num2str(size(F))]);

P = (K - omega^2*M) \ F;  % 求解得到列向量

% 列向量转矩阵（注意维度顺序）
P_matrix = reshape(P, nely+1, nelx+1); 

% 绘制彩色图
figure;
imagesc(P_matrix);        % 自动缩放颜色范围
axis equal tight;         % 等比例坐标轴+紧凑显示
axis xy;                  % 原点在左下角
colorbar;                 % 显示颜色条
colormap('jet');          % 使用jet色图（可选其他色图）
title('气压');

