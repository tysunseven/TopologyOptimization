%% 参数设置
% 几何参数
Lx = 2; Ly = Lx; nelx = 30; nely = nelx; h = Lx/nelx;
% 材料参数
rho1 = 1; rho2 = 2; E1 = 4; E2 = 20; nu = 0.34; mu1 = E1/(2*(1+nu)); mu2 = E2/(2*(1+nu));
% 优化相关参数
mode = 1; p = 48; q = -p; x = zeros(nely,nelx); % 
num_random_points = 32;

kpointinfo = struct(...
    'mu_x',     0, ...          % 波矢x分量 (标量)
    'mu_y',     0, ...           % 波矢y分量 (标量)
    'eigenvalues', zeros(num_modes, 1) ...  % 特征值数组 (num_modes × 1)
);

[boundary_mu_x, boundary_mu_y, path_distance] = generate_band_path_triangle();

boundary_points = repmat(kpointinfo, numel(boundary_mu_x), 1);
for i = 1:numel(boundary_mu_x)
    boundary_points(i).mu_x = boundary_mu_x(i); boundary_points(i).mu_y = boundary_mu_y(i);
end

random_points = repmat(kpointinfo, num_random_points, 1);

% 特征值存储
num_modes = 5; eigenvalues = zeros(length(boundary_mu_x)+num_random_points, num_modes);

%% 图形窗口
fig = figure('Position', [0 200 1400 600], 'Color', 'w', 'Visible', 'on');
% ========================= 左图初始化 =========================
subplot(1,2,1);
colors = lines(num_modes);
hold on
bandLines = gobjects(num_modes, 1); % 预创建图形对象
for l = 1:num_modes
    bandLines(l) = plot(nan, nan, 'Color', colors(l,:), 'LineWidth', 1.5);
end

% 添加带隙虚线对象 (新增代码)
gapLines = gobjects(2,1);
gapLines(1) = plot(nan, nan, '--', 'Color', [0.7 0.2 0.2], 'LineWidth', 1.5); % 下边界
gapLines(2) = plot(nan, nan, '--', 'Color', [0.7 0.2 0.2], 'LineWidth', 1.5); % 上边界


% 坐标轴设置
xticks([0, 0.292893, 0.585786, 1]) 
xticklabels({'Γ','X','M','Γ'})
xlim([0 1])
ylim([0 10])
ylabel('Frequency (Hz)')
title('Brillouin Zone Band Structure')
grid on
set(gca, 'FontSize', 12, 'FontWeight', 'bold')

% ========================= 右图初始化 =========================
subplot(1,2,2);
densityImage = imagesc(real(1-x)); % 初始化密度图
colormap(gray);
caxis([0 1]);
axis equal tight
axis off
title('Density Plot');
set(gca, 'FontSize', 12, 'FontWeight', 'bold')

%% 有限元预备
nodenrsa = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVeca = reshape(nodenrsa(1:end-1,1:end-1)+nely+1,nelx*nely,1);
edofMata = repmat(edofVeca,1,4)+repmat([0 -nely-1 -nely 1],nelx*nely,1);

Ke = [4 -1 -2 -1; -1 4 -1 -2; -2 -1 4 -1; -1 -2 -1 4]/6;
Me = h^2*[4 2 1 2; 2 4 2 1; 1 2 4 2; 2 1 2 4]/36;

iIndex = reshape(kron(edofMata,ones(4,1))',16*nelx*nely,1);
jIndex = reshape(kron(edofMata,ones(1,4))',16*nelx*nely,1);

[row, col, fixT] = init_trans(nelx,nely);

%% MMA参数
loop = 0;
change = 1;
f_history = [];
m     = 1;                % The number of general constraints. Volumn constrain
n     = nelx*nely;        % The number of design variables x_j.
xmin  = zeros( n, 1 );    % Column vector with the lower bounds for the variables x_j.
xmax  = ones( n, 1 );     % Column vector with the upper bounds for the variables x_j.
xold1 = x(:);          % xval, one iteration ago (provided that iter>1).
xold2 = x(:);          % xval, two iterations ago (provided that iter>2).
low   = [  ];             % Column vector with the lower asymptotes from the previous iteration (provided that iter>1).
upp   = [  ];             % Column vector with the upper asymptotes from the previous iteration (provided that iter>1).
a0    = 1;                % The constants a_0 in the term a_0*z.
a     = zeros(m,1);       % Column vector with the constants a_i in the terms a_i*z.
c_MMA = 1e-3*ones(m,1);   % Column vector with the constants c_i in the terms c_i*y_i.
d     = zeros(m,1);       % Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
move  = 0.2;              % Max allowed change in design variables

%% 循环体
while change>0.01 || loop<20
    loop = loop + 1; xmin = max( x(:) - move, 0 ); xmax = min( x(:) + move, 1 );

    % 局部刚度矩阵
    mu = mu1*(1-x)+mu2*x;
    rho = rho1*(1-x)+rho2*x;
    sKa = reshape(Ke(:)*(mu(:))',16*nelx*nely,1);
    sMa = reshape(Me(:)*(rho(:))',16*nelx*nely,1);
    K = sparse(iIndex,jIndex,sKa); K = (K+K')/2;
    M = sparse(iIndex,jIndex,sMa); M = (M+M')/2;
    
    random_mu_x = pi * rand(num_random_points, 1); random_mu_y = pi * rand(num_random_points, 1);
    for i = 1:numel(num_random_points)
        random_points(i).mu_x = random_mu_x(i); random_points(i).mu_y = random_mu_y(i);
    end

    combined_points = [boundary_points; random_points];

    % 并行求解能带
    dcnmax = zeros(nely, nelx);
    dcnplus1min = zeros(nely, nelx);

    parfor i = 1:length(combined_points)
        T = create_T(combined_points(i).mu_x, combined_points(i).mu_y, nelx, nely, row, col, fixT);
        K_tilde = T' * K * T; M_tilde = T' * M * T;
        [V, D] = eigs(K_tilde, M_tilde, num_modes, 'sm');
        eigenvalues(i,:) = sort(sqrt(abs(real(diag(D)))));
        combined_points(i).eigenvalues = sort(sqrt(abs(real(diag(D)))));
        phi = T*V(:,mode);  psi = T*V(:,mode+1);
        ceKmax = (mu2-mu1)*reshape(sum((phi(edofMata)*Ke).*phi(edofMata),2),nely,nelx);
        ceMmax = (rho2-rho1)*reshape(sum((phi(edofMata)*Me).*phi(edofMata),2),nely,nelx);
        dcnmax = dcnmax + sqrt(real(D(mode,mode)))^(p-1)*real((ceKmax-D(mode,mode)*ceMmax)/(2*sqrt(real(D(mode,mode)))*phi'*M*phi));
        ceKmin = (mu2-mu1)*reshape(sum((psi(edofMata)*Ke).*psi(edofMata),2),nely,nelx);
        ceMmin = (rho2-rho1)*reshape(sum((psi(edofMata)*Me).*psi(edofMata),2),nely,nelx);
        dcnplus1min = dcnplus1min + sqrt(real(D(mode+1,mode+1)))^(q-1)*real((ceKmin-D(mode+1,mode+1)*ceMmin)/(2*sqrt(real(D(mode+1,mode+1)))*psi'*M*psi));
    end
    all_eigenvalues = vertcat(combined_points.eigenvalues);
    
    boundary_points = combined_points(1:numel(boundary_mu_x));  % 提取原始边界点
    random_points = combined_points(numel(boundary_mu_x)+1:end);  % 提取随机点

    dcnmax = dcnmax * sum(eigenvalues(:, mode).^p)^(1/p-1);
    dcnplus1min = dcnplus1min * sum(eigenvalues(:, mode+1).^q).^(1/q-1);
    c = sum(eigenvalues(:, mode).^p).^(1/p)-sum(eigenvalues(:, mode+1).^q).^(1/q)+20;
    dc = dcnmax-dcnplus1min;
    % MMA优化
    
    f0val = c;                                    % scalar
    df0dx = dcnmax(:);                                % column 
    dv = ones(nely,nelx); 
    fval  = (sum(x(:))/(1000*n) - 1);            % scalar volumn constrain
    dfdx  = dv(:)'/ (0.5*n);                  % row
    [ xmma, ~, ~, ~, ~, ~, ~, ~, ~, low, upp ] = ...
        mmasub( m, n, loop, x(:), xmin, xmax, xold1, xold2, f0val,df0dx,fval,dfdx,low,upp,a0,a,c_MMA,d);
    xold2    = xold1(:);
    xold1    = x(:);
    x  = reshape(xmma,nely,nelx);
       
    change = max(abs(x(:)-xold1(:)));
      
    eigenvalues_matrix = [boundary_points.eigenvalues]; 
    % 更新图片
    for l = 1:num_modes
        set(bandLines(l), 'XData', path_distance, 'YData', eigenvalues_matrix(l, :));
    end
    % 计算带隙范围 (新增代码)
    lower_bound = max(eigenvalues(:,mode));  % 第一条能带的最大值
    upper_bound = min(eigenvalues(:,mode+1));  % 第二条能带的最小值
    
    % 更新带隙虚线 (新增代码)
    if lower_bound < upper_bound
        set(gapLines(1), 'XData', [0 1], 'YData', [lower_bound lower_bound]);
        set(gapLines(2), 'XData', [0 1], 'YData', [upper_bound upper_bound]);
    else
        set(gapLines(1), 'XData', [], 'YData', []);
        set(gapLines(2), 'XData', [], 'YData', []);
    end
    
    set(densityImage, 'CData', real(1-x));
    drawnow limitrate
   
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,upper_bound-lower_bound, ...
        mean(x(:)),change);
end

