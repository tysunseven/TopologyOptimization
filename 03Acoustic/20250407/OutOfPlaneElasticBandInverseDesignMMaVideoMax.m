%% 参数设置
% 几何参数
Lx = 2; Ly = Lx; nelx = 50; nely = nelx; h = Lx/nelx;
% 材料参数
rho1 = 1; rho2 = 2; E1 = 4; E2 = 20; nu = 0.34; mu1 = E1/(2*(1+nu)); mu2 = E2/(2*(1+nu));
% 优化相关参数
mode = 1; p = 48; q = -p; x = rand(nely, nelx); %
% 特征值存储
num_modes = 5; eigenvalues = zeros(length(mu_x_all), num_modes);
% 路径设置
x1 = linspace(0, pi, 20); y1 = zeros(1, 20); d1 = linspace(0, pi, 20);
x2 = pi * ones(1, 20); y2 = linspace(0, pi, 20); d2 = linspace(pi, 2*pi, 20);
x3 = linspace(pi, 0, 28); y3 = linspace(pi, 0, 28); d3 = linspace(2*pi, 2*pi + hypot(pi, pi), 28);
mu_x_all = [x1, x2, x3]; mu_y_all = [y1, y2, y3]; path_distance = [d1, d2, d3];
path_distance = (path_distance - min(path_distance)) / (max(path_distance) - min(path_distance));


% %% 视频设置
% videoName = 'optimization_process.mp4';
% videoObj = VideoWriter(videoName, 'MPEG-4'); % 使用MP4格式
% videoObj.FrameRate = 10;                     % 帧率10fps
% videoObj.Quality = 90;                       % 图像质量（0-100）
% open(videoObj);

%% 图形窗口
fig = figure('Position', [0 200 1400 600], 'Color', 'w');
set(fig, 'Visible', 'on'); % 设为可见以便实时观察

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
ylim([0 60])
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

Kae = [4 -1 -2 -1; -1 4 -1 -2; -2 -1 4 -1; -1 -2 -1 4]/6;
Mae = h^2*[4 2 1 2; 2 4 2 1; 1 2 4 2; 2 1 2 4]/36;

iIndexaa = reshape(kron(edofMata,ones(4,1))',16*nelx*nely,1);
jIndexaa = reshape(kron(edofMata,ones(1,4))',16*nelx*nely,1);

row = [1,((nely+1)+1):(nely+1):((nelx-1)*(nely+1)+1),nelx*(nely+1)+1,(nelx*(nely+1)+2):((nelx+1)*(nely+1))];
col = [nely+1,(2*(nely+1):(nely+1):(nelx*(nely+1))),nely+1,2:nely,nely+1];
fixdofs = setdiff(1:(nely+1)*(nelx+1),row);
fixT = sparse(fixdofs,fixdofs,ones(size(fixdofs)),(nely+1)*(nelx+1),(nely+1)*(nelx+1));

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
    loop = loop + 1;

    % 局部刚度矩阵
    mu = mu1*(1-x)+mu2*x;
    rho = rho1*(1-x)+rho2*x;
    sKa = reshape(Kae(:)*(mu(:))',16*nelx*nely,1);
    sMa = reshape(Mae(:)*(rho(:))',16*nelx*nely,1);
    Ka = sparse(iIndexaa,jIndexaa,sKa); Ka = (Ka+Ka')/2;
    Ma = sparse(iIndexaa,jIndexaa,sMa); Ma = (Ma+Ma')/2;
    
    % 并行求解能带
    dcnmax = zeros(nely, nelx);
    dcnplus1min = zeros(nely, nelx);

    parfor i = 1:length(mu_x_all)
        nu_x = mu_x_all(i);
        nu_y = mu_y_all(i);
        valcycle = [repmat(exp(-1i*nu_y),1,nelx),exp(-1i*nu_x)*exp(-1i*nu_y),repmat(exp(-1i*nu_x),1,nely)];
        Tcycle=sparse(row,col,valcycle,(nely+1)*(nelx+1),(nely+1)*(nelx+1))+fixT;
        Tcycle = Tcycle(:, ~ismember(1:size(Tcycle, 2), row));
        K_redcycle = Tcycle' * Ka * Tcycle;
        M_redcycle = Tcycle' * Ma * Tcycle;
        [Vcycle, Dcycle] = eigs(K_redcycle, M_redcycle, num_modes, 'sm');
        eigenvalues(i,:) = sort(real(diag(Dcycle)));
        phi = Tcycle*Vcycle(:,mode);  psi = Tcycle*Vcycle(:,mode+1);
        ceKmax = (mu2-mu1)*reshape(sum((phi(edofMata)*Kae).*phi(edofMata),2),nely,nelx);
        ceMmax = (rho2-rho1)*reshape(sum((phi(edofMata)*Mae).*phi(edofMata),2),nely,nelx);
        dcnmax = dcnmax + real(Dcycle(mode,mode))^(p-1)*real((ceKmax-Dcycle(mode,mode)*ceMmax)/(phi'*Ma*phi));
        ceKmin = (mu2-mu1)*reshape(sum((psi(edofMata)*Kae).*psi(edofMata),2),nely,nelx);
        ceMmin = (rho2-rho1)*reshape(sum((psi(edofMata)*Mae).*psi(edofMata),2),nely,nelx);
        dcnplus1min = dcnplus1min + real(Dcycle(mode+1,mode+1))^(q-1)*real((ceKmin-Dcycle(mode+1,mode+1)*ceMmin)/(psi'*Ma*psi));
    end

    dcnmax = dcnmax * sum(eigenvalues(:, mode).^p).^(1/p-1);
    dcnplus1min = dcnplus1min * sum(eigenvalues(:, mode+1).^q).^(1/q-1);
    c = sum(eigenvalues(:, mode).^p).^(1/p)-sum(eigenvalues(:, mode+1).^q).^(1/q)+20;
    dc = dcnmax-dcnplus1min;
    % MMA优化
    xmin = max( x(:) - move, 0 );                 
    xmax = min( x(:) + move, 1 );
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
      
    % 更新图片
    for l = 1:num_modes
        set(bandLines(l), 'XData', path_distance, 'YData', eigenvalues(:,l));
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
    % writeVideo(videoObj, getframe(fig));
end

% %% 关闭视频文件
% close(videoObj);
% disp(['视频已保存为：' videoName]);