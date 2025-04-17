%% 参数设置
Lx = 2;
Ly = 2;
nelx = 100;
nely = 100;
h = Lx/nelx;

rho1 = 1000;
E1 = 4;
nu = 0.34;
mu1 = E1/(2*(1+nu));

%% 设计变量的初始化
x = ones(nely, nelx);
rho = rho1*x;
mu = mu1*x;

%% 物理场的编号
nodenrsa = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx); % 给(nely+1)*(nelx+1)个节点编号
edofVeca = reshape(nodenrsa(1:end-1,1:end-1)+nely+1,nelx*nely,1); % 记录每个单元的右上节点的编号
edofMata = repmat(edofVeca,1,4)+repmat([0 -nely-1 -nely 1],nelx*nely,1); % 记录每个单元的四个声压自由度的编号

%% 局部矩阵的计算
Kae = [4 -1 -2 -1; -1 4 -1 -2; -2 -1 4 -1; -1 -2 -1 4]/6;
Mae = h^2*[4 2 1 2; 2 4 2 1; 1 2 4 2; 2 1 2 4]/36;

iIndexaa = reshape(kron(edofMata,ones(4,1))',16*nelx*nely,1);
jIndexaa = reshape(kron(edofMata,ones(1,4))',16*nelx*nely,1);

%% 全局矩阵的拼装
sKa = reshape(Kae(:)*(mu(:))',16*nelx*nely,1);
sMa = reshape(Mae(:)*(rho(:))',16*nelx*nely,1);
Ka = sparse(iIndexaa,jIndexaa,sKa); Ka = (Ka+Ka')/2;
Ma = sparse(iIndexaa,jIndexaa,sMa); Ma = (Ma+Ma')/2;

%% 周期性边界条件
% 完整修改后的代码
num_samples = 20;
mu_x_values = linspace(0, pi, num_samples);  % 修改点1：μx变化
mu_y_values = zeros(1, num_samples);          % 修改点2：μy固定为0
num_modes = 5;
eigenvalues = zeros(num_samples, num_modes);

for i = 1:num_samples
    mu_x = mu_x_values(i);
    mu_y = mu_y_values(i);  % 注意参数顺序
    
    T_reduced = createTransferMatrix(nelx, nely, mu_x, mu_y);
    Ka_reduced = T_reduced' * Ka * T_reduced;
    Ma_reduced = T_reduced' * Ma * T_reduced;
    
    [V, D] = eigs(Ka_reduced, Ma_reduced, num_modes, 'sm');
    omega_squared = real(diag(D));
    eigenvalues(i, :) = sort(omega_squared(1:num_modes));
end

% 绘图部分修改
figure;
hold on;
for mode = 1:num_modes
    plot(mu_x_values, sqrt(eigenvalues(:, mode))/(2*pi),...  % 修改点3：x轴数据
         'o-', 'LineWidth', 1.5);
end
xlabel('\mu_x (0 \rightarrow \pi)');        % 修改点4：坐标轴标签
title('Band Structure along \Gamma-M Path'); % 修改点5：路径名称





function T_reduced = createTransferMatrix(nelx, nely, mu_x, mu_y)
    totalNodes = (nely + 1) * (nelx + 1);
    rows = [];
    cols = [];
    vals = [];
    
    % 处理上边界节点（左端和中间nelx-1个）
    % 上边左端节点
    rows = [1];
    cols = [nely + 1];
    vals = [exp(-1i*mu_y)];
    
    % 上边中间nelx-1个节点
    for i = 1:(nelx-1)
        row_i = 1 + i*(nely+1);
        col_i = (i+1)*(nely+1);
        rows = [rows, row_i];
        cols = [cols, col_i];
        vals = [vals, exp(-1i*mu_y)];
    end
    
    % 处理右上角节点
    row_ur = nelx*(nely+1) + 1;
    col_ur = nely + 1;
    rows = [rows, row_ur];
    cols = [cols, col_ur];
    vals = [vals, exp(-1i*(mu_x + mu_y))];
    
    % 处理右侧中间nely-1个节点
    if nely > 1
        rows_right = (nelx*(nely+1)+2) : (nelx*(nely+1)+nely);
        cols_right = 2:nely;
        rows = [rows, rows_right];
        cols = [cols, cols_right];
        vals = [vals, exp(-1i*mu_x)*ones(1,nely-1)];
    end
    
    % 处理右下角节点
    row_br = (nelx+1)*(nely+1);
    col_br = nely + 1;
    rows = [rows, row_br];
    cols = [cols, col_br];
    vals = [vals, exp(-1i*(mu_x + mu_y))];
    
    % 构建稀疏矩阵
    T = sparse(rows, cols, vals, totalNodes, totalNodes);
    
    % 添加对角线元素（未被处理的节点）
    allNodes = 1:totalNodes;
    handledNodes = unique(rows);
    unhandledNodes = setdiff(allNodes, handledNodes);
    T = T + sparse(unhandledNodes, unhandledNodes, ones(length(unhandledNodes),1), totalNodes, totalNodes);
    
    % 删除全零列
    zeroCols = all(T == 0, 1);
    T_reduced = T(:, ~zeroCols);
end