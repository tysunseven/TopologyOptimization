%% 参数设置
%x = zeros(nely, nelx);
% x1 =  rand(10, 10);
% x2 = zeros(10, 10);
% x = [x1,x2; x2,x1']
x = readmatrix('micro2.xlsx');
figure('Position', [100 100 400 400]);
colormap(gray);
imagesc(1 - x);
clim([0 1]);
axis equal tight;
axis off;
title('Material Distribution');
drawnow;

Lx = 2;
Ly = 2;
nelx = 20;
nely = 20;
h = Lx/nelx;

rho1 = 1;
rho2 = 2;
E1 = 4;
E2 = 20;
nu = 0.34;
mu1 = E1/(2*(1+nu));
mu2 = E2/(2*(1+nu));

%% 有限元

rho = rho1*(1-x)+rho2*x;
mu = mu1*(1-x)+mu2*x;

nodenrsa = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVeca = reshape(nodenrsa(1:end-1,1:end-1)+nely+1,nelx*nely,1);
edofMata = repmat(edofVeca,1,4)+repmat([0 -nely-1 -nely 1],nelx*nely,1);

Kae = [4 -1 -2 -1; -1 4 -1 -2; -2 -1 4 -1; -1 -2 -1 4]/6;
Mae = h^2*[4 2 1 2; 2 4 2 1; 1 2 4 2; 2 1 2 4]/36;

iIndexaa = reshape(kron(edofMata,ones(4,1))',16*nelx*nely,1);
jIndexaa = reshape(kron(edofMata,ones(1,4))',16*nelx*nely,1);

sKa = reshape(Kae(:)*(mu(:))',16*nelx*nely,1);
sMa = reshape(Mae(:)*(rho(:))',16*nelx*nely,1);

Ka = sparse(iIndexaa,jIndexaa,sKa); Ka = (Ka+Ka')/2;
Ma = sparse(iIndexaa,jIndexaa,sMa); Ma = (Ma+Ma')/2;

%% 能带结构计算 (新增部分)
[mu_x_all, mu_y_all, path_distance] = generate_band_path_triangle();

% 特征值存储
num_modes = 5;
eigenvalues = zeros(length(mu_x_all), num_modes);

row = [1,((nely+1)+1):(nely+1):((nelx-1)*(nely+1)+1),nelx*(nely+1)+1,(nelx*(nely+1)+2):((nelx+1)*(nely+1))];
col = [nely+1,(2*(nely+1):(nely+1):(nelx*(nely+1))),nely+1,2:nely,nely+1];
fixdofs = setdiff(1:(nely+1)*(nelx+1),row);
fixT = sparse(fixdofs,fixdofs,ones(size(fixdofs)),(nely+1)*(nelx+1),(nely+1)*(nelx+1));

% 主计算循环
parfor i = 1:length(mu_x_all)
    nu_x = mu_x_all(i);
    nu_y = mu_y_all(i);
    val = [repmat(exp(-1i*nu_y),1,nelx),exp(-1i*nu_x)*exp(-1i*nu_y),repmat(exp(-1i*nu_x),1,nely)];
    T=sparse(row,col,val,(nely+1)*(nelx+1),(nely+1)*(nelx+1))+fixT;
    T = T(:, ~ismember(1:size(T, 2), row));
    K_red = T' * Ka * T;
    M_red = T' * Ma * T;
    [~, D] = eigs(K_red, M_red, num_modes, 'sm');
    eigenvalues(i,:) = sort(real(diag(D)));
end

p = 48;                                 % 示例 p 值，根据实际需求调整
f = sum(eigenvalues(:, 1).^p).^(1/p);   % 手动计算 p-范数
disp(max(eigenvalues(:, 1)));
disp(f);

% mu_x = pi*0.9;
% mu_y = pi*0.8;
% val = [repmat(exp(-1i*mu_y),1,nelx),exp(-1i*mu_x)*exp(-1i*mu_y),repmat(exp(-1i*mu_x),1,nely)];
% T=sparse(row,col,val,(nely+1)*(nelx+1),(nely+1)*(nelx+1))+fixT;
% T = T(:, ~ismember(1:size(T, 2), row));
% K_red = T' * Ka * T;
% M_red = T' * Ma * T;
% [V, D] = eigs(K_red, M_red, 5, 'sm');  % 求最小特征值及对应特征向量
% disp(D);

%% 能带结构可视化（原plot_band_structure功能）
figure('Position', [500 100 800 600]);
colors = lines(size(eigenvalues,2));

hold on
for m = 1:size(eigenvalues,2)
    plot(path_distance, eigenvalues(:,m),...
        'Color', colors(m,:), 'LineWidth', 1.5)
end

% 原标注 (需修改为) 正方形
% xticks([0, 0.25, 0.5, 0.75 1]) 
% xticklabels({'Γ','X','M','','Γ'})

%三角形
xticks([0, 0.292893, 0.585786, 1]) 
xticklabels({'Γ','X','M','Γ'})


xlim([0 1])

ylabel('Frequency (Hz)')
title('Brillouin Zone Band Structure')
grid on
set(gca, 'FontSize', 12, 'FontWeight', 'bold')