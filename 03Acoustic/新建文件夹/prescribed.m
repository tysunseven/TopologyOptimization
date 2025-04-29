%% 参数设置
Lx = 2;
Ly = 2;
nelx = 20;
nely = 20;
h = Lx/nelx;
mu_x = pi;
mu_y = pi;
volfrac=0.4;

rho1 = 1;
rho2 = 2;
E1 = 4;
E2 = 20;
nu = 0.34;
mu1 = E1/(2*(1+nu));
mu2 = E2/(2*(1+nu));

x = readmatrix('micro2.xlsx');

omegac=9.5;
omegad=0.5;

%% 能带结构计算 (新增部分)
% 第一段 Γ→X (0,0)->(pi,0)
x1 = linspace(0, pi, 20);
y1 = zeros(1, 20);
d1 = linspace(0, pi, 20);

% 第二段 X→M (pi,0)->(pi,pi)
x2 = pi * ones(1, 20);
y2 = linspace(0, pi, 20);
d2 = linspace(pi, 2*pi, 20);

% 第三段 M→Γ (pi,pi)->(0,0) 直线
x3 = linspace(pi, 0, 28);
y3 = linspace(pi, 0, 28);
seg3_length = hypot(pi, pi); % 直线距离为pi*sqrt(2)
d3 = linspace(2*pi, 2*pi + seg3_length, 28);

% 合并所有数据
mu_x_all = [x1, x2, x3];
mu_y_all = [y1, y2, y3];
path_distance = [d1, d2, d3];

% 归一化路径参数
path_distance = (path_distance - min(path_distance)) / ...
               (max(path_distance) - min(path_distance));

% 特征值存储
num_modes = 5;
eigenvalues = zeros(length(mu_x_all), num_modes);

%% 有限元
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
c_MMA = 1000*ones(m,1);   % Column vector with the constants c_i in the terms c_i*y_i.
d     = zeros(m,1);       % Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
move  = 0.2;              % Max allowed change in design variables
while change>0.01
    loop = loop + 1;

    mu = mu1*(1-x)+mu2*x;
    rho = rho1*(1-x)+rho2*x;
    sKa = reshape(Kae(:)*(mu(:))',16*nelx*nely,1);
    sMa = reshape(Mae(:)*(rho(:))',16*nelx*nely,1);
    Ka = sparse(iIndexaa,jIndexaa,sKa); Ka = (Ka+Ka')/2;
    Ma = sparse(iIndexaa,jIndexaa,sMa); Ma = (Ma+Ma')/2;

    val = [repmat(exp(-1i*mu_y),1,nelx),exp(-1i*mu_x)*exp(-1i*mu_y),repmat(exp(-1i*mu_x),1,nely)];
    T=sparse(row,col,val,(nely+1)*(nelx+1),(nely+1)*(nelx+1))+fixT;
    T = T(:, ~ismember(1:size(T, 2), row));
    K_red = T' * Ka * T; M_red = T' * Ma * T;
    [V, D] = eigs(K_red, M_red, 1, 'sm');  % 求最小特征值及对应特征向量
    f = real(D(1,1));                            % 特征值存入 f（标量）
    omega_prescribed = real(D(1,1));
    disp(omega_prescribed);
    phi = T*V(:,1);                          % 特征向量存入 phi（列向量）


    ceK = (mu2-mu1)*reshape(sum((phi(edofMata)*Kae).*phi(edofMata),2),nely,nelx);
    ceM = (rho2-rho1)*reshape(sum((phi(edofMata)*Mae).*phi(edofMata),2),nely,nelx);
    dc = real((ceK-f*ceM)/(phi'*Ma*phi));

    f_history = [f_history f];  % Store objective function value


    xmin = max( x(:) - move, 0 );                 
    xmax = min( x(:) + move, 1 );
    f_prescribed = exp(-(((omega_prescribed-omegac)/omegad)^2)/2);
    disp(f_prescribed);
    f0val = f_prescribed;                                    % scalar
    df0dx = -f_prescribed*((omega_prescribed-omegac)/omegad^2)*dc(:);                                % column 
    dv = ones(nely,nelx); 
    fval  = -(sum(x(:))/(volfrac*n) - 1);            % scalar volumn constrain
    dfdx  = -dv(:)'/ (volfrac*n);                  % row
    [ xmma, ~, ~, ~, ~, ~, ~, ~, ~, low, upp ] = ...
        mmasub( m, n, loop, x(:), xmin, xmax, xold1, xold2, f0val,df0dx,fval,dfdx,low,upp,a0,a,c_MMA,d);
    % Update MMA Variables
    xold2    = xold1(:);
    xold1    = x(:);
    x   = reshape(xmma,nely,nelx);
       
    change = max(abs(x(:)-xold1(:)));
      
    %% PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,f, ...
        mean(x(:)),change);

    parfor i = 1:length(mu_x_all)
        nu_x = mu_x_all(i);
        nu_y = mu_y_all(i);
        valcycle = [repmat(exp(-1i*nu_y),1,nelx),exp(-1i*nu_x)*exp(-1i*nu_y),repmat(exp(-1i*nu_x),1,nely)];
        Tcycle=sparse(row,col,valcycle,(nely+1)*(nelx+1),(nely+1)*(nelx+1))+fixT;
        Tcycle = Tcycle(:, ~ismember(1:size(Tcycle, 2), row));
        K_redcycle = Tcycle' * Ka * Tcycle;
        M_redcycle = Tcycle' * Ma * Tcycle;
        [~, Dcycle] = eigs(K_redcycle, M_redcycle, num_modes, 'sm');
        eigenvalues(i,:) = sort(real(diag(Dcycle)));
    end
    
    figure('Position', [0 100 1400 600]); % 加宽窗口适应两张图

% ========================= 左图：能带结构 =========================
subplot(1,2,1);
colors = lines(size(eigenvalues,2));

hold on
for l = 1:size(eigenvalues,2)
    plot(path_distance, eigenvalues(:,l),...
        'Color', colors(l,:), 'LineWidth', 1.5)
end

% 设置坐标轴
xticks([0, 0.292893, 0.585786, 1]) 
xticklabels({'Γ','X','M','Γ'})
xlim([0 1])
ylim([0 60])

% 标签样式
ylabel('Frequency (Hz)')
title('Brillouin Zone Band Structure')
grid on
set(gca, 'FontSize', 12, 'FontWeight', 'bold')

% ========================= 右图：密度图 =========================
subplot(1,2,2);
colormap(gray); 
imagesc(real(1-x)); 
caxis([0 1]); 
axis equal; 
axis off; 
title('Density Plot');
set(gca, 'FontSize', 12, 'FontWeight', 'bold')

drawnow;
    %filename = sprintf('loop_%d.xlsx', loop);
    %writematrix(x, filename);
end