% 参数设置
nelx = 100;  % x方向网格数
nely = 50;  % y方向网格数
rho1 = 1.225;  % 空气的密度1.225kg/m^3
kappa1 = 141834.999;  % 空气的体积模量
xi = ones(nely, nelx);

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

iIndex = reshape(kron(edofMat,ones(4,1))',16*nelx*nely,1);
jIndex = reshape(kron(edofMat,ones(1,4))',16*nelx*nely,1);

sK = reshape(KE(:)*xi(:)',16*nelx*nely,1);
sM = reshape(ME(:)*xi(:)',16*nelx*nely,1);

K = sparse(iIndex,jIndex,sK); K = (K+K')/2;
M = sparse(iIndex,jIndex,sM); M = (M+M')/2;

K = K/rho1;
M = M/kappa1;

% 求解前 10 个最小特征值（声学问题通常关注低频模态）
num_modes = 10;
[V, D] = eigs(K, M, num_modes, 'sm'); 

% 提取特征值 omega^2（按升序排列）
omega_squared = diag(D);

f = sqrt(omega_squared) / (2*pi);

disp(f);

Lx = 2;
Ly = 1;
c = 340.27;

max_mode = 10; % 生成m和n从1到10的组合

% 创建所有可能的(m,n)对
[m_grid, n_grid] = meshgrid(0:max_mode, 0:max_mode);
m = m_grid(:);
n = n_grid(:);

% 计算理论频率
frequencies = (c / 2) * sqrt( (m./Lx).^2 + (n./Ly).^2 );

% 排序频率
sorted_frequencies = sort(frequencies);

% 取前10个
top10 = sorted_frequencies(1:10)/50;

% 显示结果
disp('前10个理论特征频率 (Hz):');
disp(top10);

