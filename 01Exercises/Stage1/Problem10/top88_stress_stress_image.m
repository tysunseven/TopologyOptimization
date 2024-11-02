%% Set Parameters

nelx=40;
nely=20;
volfrac=0.3;
penal=3.0;
rmin=nely*0.06;

ft=3;

%% MATERIAL PROPERTIES

E0 = 100;
Emin = 1e-9*E0;
nu = 0.3;
In_spring=0.1;
Out_spring=0.1;
a = 0.5;
b = 0.5;
pnorm = 12;

B = [1/2 0 -1/2 0 -1/2 0 1/2 0;
     0 1/2 0 1/2 0 -1/2 0 -1/2;
     1/2 1/2 1/2 -1/2 -1/2 -1/2 -1/2 1/2];

D = (E0 / (1 - nu^2)) * [1, nu, 0;
                        nu, 1, 0;
                        0, 0, (1 - nu) / 2];

V = [1, -1/2, 0;
    -1/2, 1, 0;
    0, 0, 3];

Q = B'*D'*V*D*B;

%% PREPARE FINITE ELEMENT ANALYSIS

A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);

edof = zeros(nelx * nely, 8);
elem_num = 0;

for j = 1:nelx
    for i = 1:nely
        elem_num = elem_num + 1;

        n1 = (j) * (nely + 1) + i;
        n2 = (j - 1) * (nely + 1) + i;
        n3 = (j - 1) * (nely + 1) + i + 1;
        n4 = (j) * (nely + 1) + i + 1;

        edof(elem_num, :) = [2 * n1 - 1, 2 * n1, 2 * n2 - 1, 2 * n2, ...
                             2 * n3 - 1, 2 * n3, 2 * n4 - 1, 2 * n4];
    end
end


%% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)

U = zeros(2*(nely+1)*(nelx+1),1);
Lambda = zeros(2*(nely+1)*(nelx+1),1); % compliant add to this
L= zeros(2*(nely+1)*(nelx+1),1);       % compliant add to this

% BC with symmetry 上半部分
%din = (nely+1)*2-1; dout = 2*(nely+1)*(nelx+1)-1;
%fixeddofs =union([1,2],2*(nely+1):2*(nely+1):2*(nely+1)*(nelx+1));%fix x of top left and y of bottom edge

% BC with symmetry 下半部分
din = 1; dout = 2*(nely+1)*nelx+1;
fixeddofs =union([2*(nely+1)-1,2*(nely+1)],2:2*(nely+1):2*(nely+1)*nelx+2);%fix x of top left and y of bottom edge


% BC without symmetry
%din=2*ceil((nely+1)/2)-1;dout = 2*(nelx)*(nely+1)+2*ceil((nely+1)/2)-1;
%fixeddofs = union([1:2], [2*(nely+1)-1,2*(nely+1)]);

F = sparse(din,1,1,2*(nely+1)*(nelx+1),1);
alldofs = 1:2*(nely+1)*(nelx+1);
freedofs = setdiff(alldofs,fixeddofs);
%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
    for j1 = 1:nely
        e1 = (i1-1)*nely+j1;
        for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
            for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                e2 = (i2-1)*nely+j2;
                k = k+1;
                iH(k) = e1;
                jH(k) = e2;
                sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% INITIALIZE ITERATION
xval = repmat(volfrac,nely,nelx);
vmstress = zeros(nely, nelx);
xPhys = xval;
beta = 1; % 初始Heaviside参数
xTilde = xval; % 用于Heaviside的中间变量
if ft == 1
    xPhys = xval;
elseif ft == 2
    xPhys( : ) = reshape( ( H * xval( : ) ) ./ Hs, nely, nelx );
elseif ft == 3 % Heaviside正则化
    xPhys = 1 - exp(-beta * xTilde) + xTilde * exp(-beta);
end
loop = 0;
normalized_coefficient = 1;
vmstress_pnorm = 1;
change = 1;
loopbeta = 0;
c_history = [];
max_vmstress_history = [];
normalized_vmstress_pnorm_history = [];
vmstress_pnorm_history = [];
%% INITIALIZE MMA OPTIMIZER
%  DESIGN UPDATE BY THE MMA METHOD
m     = 1;                % The number of general constraints. Volumn constrain
n     = nelx*nely;        % The number of design variables x_j.
xmin  = zeros(n,1);       % Column vector with the lower bounds for the variables x_j.
xmax  = ones(n,1);        % Column vector with the upper bounds for the variables x_j.
xold1 = xval(:);          % xval, one iteration ago (provided that iter>1).
xold2 = xval(:);          % xval, two iterations ago (provided that iter>2).
low   = [  ];             % Column vector with the lower asymptotes from the previous iteration (provided that iter>1).
upp   = [  ];             % Column vector with the upper asymptotes from the previous iteration (provided that iter>1).
a0    = 1;                % The constants a_0 in the term a_0*z.
a     = zeros(m,1);       % Column vector with the constants a_i in the terms a_i*z.
c_MMA = 1000*ones(m,1);   % Column vector with the constants c_i in the terms c_i*y_i.
d     = zeros(m,1);       % Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
move  = 0.1;              % Max allowed change in design variables

figure;

%% START ITERATION

tic;
while change > 0.003
    loopbeta = loopbeta+1;
    loop = loop + 1;
    if loop > 1
        normalized_coefficient = 0.5 * max(vmstress(:))/vmstress_pnorm+0.5*normalized_coefficient;
    end
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    % Adding external springs
    K(din,din) = K(din,din) + In_spring;
    K(dout,dout) = K(dout,dout) + Out_spring;
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);
    L(dout)=1;
    Lambda(freedofs) = -K(freedofs,freedofs)\L(freedofs);
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = reshape(sum((Lambda(edofMat)*KE).*U(edofMat),2),nely,nelx);       % ce=Lambda'*K*U
    c = U(dout);
    dc = penal*(E0-Emin)*xPhys.^(penal-1).*ce;
    dv = ones(nely,nelx);                                                  % eq 6 in top88
    c_history = [c_history c];

    Ue = U(edof);
    UeQ = Ue * Q;
argument = sum(UeQ .* Ue, 2);

% 确保非负
argument(argument < 0) = 0;

% 计算 vmstress
vmstress_vector = xPhys(:).^2 .* sqrt(argument);

    %vmstress_vector = xPhys(:).^2 .* sqrt(sum(UeQ .* Ue, 2));
    vmstress = reshape(vmstress_vector, nely, nelx);
    vmstress = double(vmstress);
    vmstress_p = vmstress.^pnorm;
    sum_vmstress_p = sum(vmstress_p(:));
    vmstress_pnorm = sum_vmstress_p^(1/pnorm);
    normalized_vmstress_pnorm = normalized_coefficient * vmstress_pnorm;
    vmstress_pnorm_history = [vmstress_pnorm_history vmstress_pnorm];
    max_vmstress_history = [max_vmstress_history max(vmstress(:))];
normalized_vmstress_pnorm_history = [normalized_vmstress_pnorm_history normalized_vmstress_pnorm];


    %% FILTERING/MODIFICATION OF SENSITIVITIES
    if ft == 1
        dc(:) = H*(xval(:).*dc(:))./Hs./max(1e-3,xval(:));
    elseif ft == 2
        dc(:) = H*(dc(:)./Hs);
        dv(:) = H*(dv(:)./Hs);
    elseif ft == 3
        dx = beta * exp(-beta * xTilde) + exp(-beta);
        dc(:) = H*(dc(:).*dx(:)./Hs);
        dv(:) = H*(dv(:).*dx(:)./Hs);
    end

    %% METHOD OF MOVING ASYMPTOTES
    xmin = max( xval(:) - move, 0 );
    xmax = min( xval(:) + move, 1 );
    f0val = c;                                    % scalar
    df0dx = dc(:);                                % column
    fval  = sum(xPhys(:))/(volfrac*n) - 1;        % scalar volumn constrain
    dfdx  = dv(:)'/ (volfrac*n);                  % row
    [ xmma, ~, ~, ~, ~, ~, ~, ~, ~, low, upp ] = ...
        mmasub( m, n, loop, xval(:), xmin, xmax, xold1, xold2, f0val,df0dx,fval,dfdx,low,upp,a0,a,c_MMA,d);
    % Update MMA Variables
    xold2    = xold1(:);
    xold1    = xval(:);
    xval     = reshape(xmma,nely,nelx);

    % FILTERING
    if ft == 1
        xPhys = xval;
    elseif ft == 2
        xPhys(:) = (H *xval(:))./ Hs;
    elseif ft == 3
        xTilde(:) = (H * xval(:))./ Hs;
        xPhys = 1 - exp(-beta * xTilde) + xTilde * exp(-beta);
    end

    if ft == 3 && beta < 512 && (loopbeta >= 50 || change <= 0.01)
        beta = 2*beta;
        loopbeta = 0;
        change = 1;
        fprintf('Parameter beta increased to %g.\n',beta);
    end

    change = max(abs(xval(:)-xold1(:)));

    
    
% for i = 1:nely
%     for j = 1:nelx
%         u_i = [U(2 * ((j+1 - 1) * (nely + 1) + i)-1),U(2 * ((j+1 - 1) * (nely + 1) + i)),U(2 * ((j - 1) * (nely + 1) + i)-1),U(2 * ((j - 1) * (nely + 1) + i)),U(2 * ((j - 1) * (nely + 1) + i+1)-1),U(2 * ((j - 1) * (nely + 1) + i+1)),U(2 * ((j+1 - 1) * (nely + 1) + i+1)-1),U(2 * ((j+1 - 1) * (nely + 1) + i+1))];  % 需要根据 U 提取出第 (i, j) 单元的位移向量
%         result = xPhys(i,j)^2*sqrt(u_i * Q * u_i');
%         vmstress2(i, j) = result;
%     end
% end

    %% PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, ...
        mean(xPhys(:)),change);
    fprintf('The value of normalized_vmstress_pnorm is: %.4f, and the maximum value in vmstress is: %.4f\n', normalized_vmstress_pnorm, max(vmstress(:)));

    
    %colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
    subplot(1, 2, 1);  % 创建 1 行 2 列的子图，并选择第 1 个位置
    colormap(gray);
    imagesc(1-xPhys);  % 绘制矩阵图像
    title('1 - xPhys');  % 添加标题
    axis equal;  % 设置坐标轴比例相同
    axis off;

    subplot(1, 2, 2);  % 选择第 2 个子图位置
    colormap(gray);
    imagesc(-vmstress);  % 绘制应力矩阵图像
    title('vmstress');  % 添加标题
    axis equal;  % 设置坐标轴比例相同
    axis off;

    

    drawnow;
end
t1 = toc; % 记录第一个时间点
disp(['代码块1的运行时间: ', num2str(t1), ' 秒']);


% 假设 max_vmstress_history、normalized_vmstress_pnorm_history 和 vmstress_pnorm_history 已经在循环中被记录

figure('Position', [100, 100, 800, 600]); % 创建一个更大的图形窗口，设置尺寸为800x600
plot(max_vmstress_history, '-o', 'DisplayName', 'Max', 'LineWidth', 0.5, 'MarkerSize', 4); % 设置线宽和标记大小
hold on;
plot(vmstress_pnorm_history, '-^', 'DisplayName', 'Pnorm', 'LineWidth', 0.5, 'MarkerSize', 4); % 设置线宽和标记大小
plot(normalized_vmstress_pnorm_history, '-s', 'DisplayName', 'Normalized Pnorm', 'LineWidth', 0.5, 'MarkerSize', 4); % 设置线宽和标记大小
hold off;

% 设置图形属性
title('Max, Pnorm and Normalized Pnorm History'); % 设置标题
xlabel('Iteration'); % 设置X轴标签
ylabel('Values'); % 设置Y轴标签
legend('Location', 'best'); % 将图例放置在最佳位置
grid on; % 显示网格

% 创建子图 2: 绘制 normalized_vmstress_pnorm_history 与 max_vmstress_history 的比值
% subplot(2,1,2); % 激活第二个子图
% ratio = normalized_vmstress_pnorm_history ./ max_vmstress_history; % 计算比值
% plot(ratio, '-d', 'DisplayName', 'Normalized Pnorm / Max VMStress'); % 绘制比值曲线，使用菱形标记
% title('Ratio of Normalized Pnorm to Max VMStress'); % 设置子图 2 的标题
% xlabel('Iteration'); % 设置X轴标签
% ylabel('Ratio'); % 设置Y轴标签
% legend show; % 显示图例
% grid on; % 显示网格

    drawnow;

save_filename = sprintf('Stress_constraints_without_constraints_iter%d_obj%.4f.png', loop, c);
saveas(gcf, save_filename);

for i = 1:nely
        for j = 1:nelx
            u_i = [U(2 * ((j+1 - 1) * (nely + 1) + i)-1),U(2 * ((j+1 - 1) * (nely + 1) + i)),U(2 * ((j - 1) * (nely + 1) + i)-1),U(2 * ((j - 1) * (nely + 1) + i)),U(2 * ((j - 1) * (nely + 1) + i+1)-1),U(2 * ((j - 1) * (nely + 1) + i+1)),U(2 * ((j+1 - 1) * (nely + 1) + i+1)-1),U(2 * ((j+1 - 1) * (nely + 1) + i+1))];  % 需要根据 U 提取出第 (i, j) 单元的位移向量
            result = xPhys(i,j)^2*sqrt(u_i * Q * u_i');
            vmstress(i, j) = result;
        end
end


% subplot(1, 2, 1);  % 创建 1 行 2 列的子图，并选择第 1 个位置
%     colormap(gray);
%     imagesc(1-xPhys);  % 绘制矩阵图像
%     title('1 - xPhys');  % 添加标题
%     axis equal;  % 设置坐标轴比例相同
%     axis off;
% 
%     % 3. 绘制第二个子图：vmstress
%     subplot(1, 2, 2);  % 选择第 2 个子图位置
%     imagesc(-vmstress);  % 绘制应力矩阵图像
%     title('vmstress');  % 添加标题
%     axis equal;  % 设置坐标轴比例相同
%     axis off;










% figure('Position', [100, 100, 1200, 500]);  % 增加图像的宽度，将整张图像拉长
% 
% % 调整结构图像的大小和位置
% subplot('Position', [0.05, 0.15, 0.4, 0.7]);  % 调整位置和大小 [左, 下, 宽, 高]
% colormap(gray);
% imagesc(1-xPhys); 
% caxis([0 1]); 
% axis equal; 
% axis off; 
% title('Optimized Structure');
% 
% % 调整目标函数历史图像的大小和位置
% subplot('Position', [0.55, 0.15, 0.4, 0.7]);  % 调整位置和大小 [左, 下, 宽, 高]
% plot(1:loop, c_history, '-o', 'LineWidth', 2);  % 使用 -o 添加数据点标记
% xlabel('Iteration');
% ylabel('Objective Function Value');
% title('Objective Function History');
% 
% % 保存图像
% save_filename = sprintf('Stress_constraints_iter%d_obj%.4f.png', loop, c);
% saveas(gcf, save_filename);

% for i = 1:nely
%         for j = 1:nelx
%             u_i = [U(compute_index(j+1,i)-1),U(compute_index(j+1,i)),U(compute_index(j,i)-1),U(compute_index(j,i)),U(compute_index(j,i+1)-1),U(compute_index(j,i+1)),U(compute_index(j+1,i+1)-1),U(compute_index(j+1,i+1))];  % 需要根据 U 提取出第 (i, j) 单元的位移向量
%             result = xPhys(i,j)^2*sqrt(u_i * Q * u_i');
%             vmstress(i, j) = result;
%         end
% end



function idx = compute_index(index_x, index_y)
    idx = 2 * ((index_x - 1) * (30 + 1) + index_y);
end


