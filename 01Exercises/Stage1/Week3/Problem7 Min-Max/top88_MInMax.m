%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE Nov, 2010 %%%%
%%%% For Min-Max formulation, multiple input for MMA%%%%
clear;close all;clc;
nelx=120;
nely=40;
volfrac=0.5;
penal=3.0;  % typically 3
rmin=3.6;  
ft=3;       
%% plotting Design and convergence
close all
%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-8;
nu = 0.3;
ndofs = 2*(nely+1)*(nelx+1);    % Number of DOF in the system
nelements = nelx*nely;          % Number of Elements in the system
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
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
total_dofs = 2*(nely+1)*(nelx+1);  % 自由度总数
multiple_load_cases = 3;           % 载荷工况数量

% 创建稀疏矩阵F
F = sparse(total_dofs, multiple_load_cases);  % 初始化大小为 (自由度数 x 载荷工况数) 的稀疏矩阵

% 定义第一列 (第一载荷工况)
F(2*(nely+1)*(nelx/4) + 2, 1) = -1; 

% 定义第二列 (第二载荷工况)
F(2*3*(nely+1)*(nelx/4)+2, 2) = -1;

% 定义第三列 (第三载荷工况)
F(2*(nely+1)*(nelx/2) + 2, 3) = -1;  


fixeddofs = [2*(nely+1)-1,2*(nely+1),2*(nelx+1)*(nely+1)-1, 2*(nelx+1)*(nely+1)]; % 左下角右下角xy分量都固定
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);
% Deformation matrix is defines (only necesarry for multiple load cases)
U = zeros(2*(nely+1)*(nelx+1),multiple_load_cases);  
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
xPhys = xval;
beta = 1;
if ft == 1                                                                 
    xPhys = xval;
elseif ft == 2
    xPhys( : ) = reshape( ( H * xval( : ) ) ./ Hs, nely, nelx );
elseif ft == 3
  xTilde = xval;
  xPhys = 1-exp(-beta*xTilde)+xTilde*exp(-beta);
end
loop = 0;
change = 1;
loopbeta = 0;
c_history = [];
%% INITIALIZE MMA OPTIMIZER
% DESIGN UPDATE BY THE MMA METHOD
% Modifed for Min Max formulation with multiple constrains
p=multiple_load_cases;
q=1;
m = 2*p+q;                % The number of general constraints. Volumn constrain
n = nelx*nely;            % The number of design variables x_j.
xmin  = zeros(n,1);       % Column vector with the lower bounds for the variables x_j.
xmax  = ones(n,1);        % Column vector with the upper bounds for the variables x_j.
xold1 = xval(:);          % xval, one iteration ago (provided that iter>1).
xold2 = xval(:);          % xval, two iterations ago (provided that iter>2).
low   = [  ];             % Column vector with the lower asymptotes from the previous iteration (provided that iter>1).
upp   = [  ];             % Column vector with the upper asymptotes from the previous iteration (provided that iter>1).
a0    = 1;                % The constants a_0 in the term a_0*z.
a     = 1*ones(m,1);      % Column vector with the constants a_i in the terms a_i*z.
a(end)=0;
c_MMA = 1000*ones(m,1);   % Column vector with the constants c_i in the terms c_i*y_i.
d     = zeros(m,1);       % Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
move  = 0.1;              % Max allowed change in design variables
fval = zeros(m,1);
dfdx = zeros(m,n);
%% START ITERATION
while change > 0.002
    loopbeta = loopbeta+1;
    loop = loop + 1;
    %% FE-ANALYSIS
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
    K = sparse(iK,jK,sK);
    K = (K+K')/2;
    U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:); 
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    % Calculating the compliance + Implement multiple load cases
    c = zeros(multiple_load_cases,1);
    dc = zeros(nely,nelx,multiple_load_cases);
    for i = 1:size(F,2) %size(,2) fucntion returns how many columns in F
        Ui = U(:,i);
        ce = reshape(sum((Ui(edofMat)*KE).*Ui(edofMat),2),nely,nelx);
        c(i) = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
        dc(:,:,i) = - penal*(E0-Emin)*xPhys.^(penal-1).*ce;
    end
    dv = ones(nely,nelx);
    c_history = [c_history c];  % Store objective function value
    z=abs(max(c));    
    %% FILTERING/MODIFICATION OF SENSITIVITIES
    %  ft=1 is sensitivity filter, ft=2 is density filter
    if ft == 1
        dc(:) = H*(xval(:).*dc(:))./Hs./max(1e-3,xval(:));
    elseif ft == 2
        for i=1:multiple_load_cases
        dc_temp = dc(:,:,i);
        dc_temp(:) = H*(dc_temp(:)./Hs);
        dc(:,:,i) = dc_temp;
        end
        dv(:) = H*(dv(:)./Hs);
        elseif ft == 3
    dx = beta*exp(-beta*xTilde)+exp(-beta);
    for i = 1:multiple_load_cases
    dc_temp = dc(:,:,i);
    dc_temp(:) = H*(dc_temp(:).*dx(:)./Hs);
    dc(:,:,i) = dc_temp;
end
    dv(:) = H*(dv(:).*dx(:)./Hs);
    end
    %% METHOD OF MOVING ASYMPTOTES
    xmin = max( xval(:) - move, 0 ); 
    xmax = min( xval(:) + move, 1 );
    % Modifed for Min Max formulation
    f0val = 0;  
    df0dx = zeros(n, 1);
    fval(1:p) = c(1:p);
    fval(p+1:2*p) = -c(1:p);
    fval(m) = sum(xPhys(:))/(volfrac*n) - 1;
    for i=1:multiple_load_cases
        dfdx(i,:) = reshape(dc(:,:,i),1,n);
        dfdx(p+i,:) = reshape(-dc(:,:,i),1,n);
    end
    dfdx(m,:) = transpose(dv(:))/(n*volfrac);
    [ xmma, ~, ~, ~, ~, ~, ~, ~, ~, low, upp ] = ...
        mmasub( m, n, loop, xval(:), xmin, xmax, xold1, xold2, f0val,df0dx,fval,dfdx,low,upp,a0,a,c_MMA,d);
    % Update MMA Variables
    xold2 = xold1(:);
    xold1 = xval(:);
    xval =  reshape(xmma,nely,nelx);
    % FILTERING
    if ft == 1         % sensitivity filter
        xPhys = xval;
    elseif ft == 2     % density filter
        xPhys(:) = (H *xval(:))./ Hs;
        elseif ft == 3
    xTilde(:) = (H * xval(:))./ Hs;
    xPhys = 1 - exp(-beta * xTilde) + xTilde * exp(-beta); % Heaviside正则化
    end        
    change = max(abs(xval(:)-xold1(:)));
    
    %% PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,z, ...
        mean(xval(:)),change);
    
    colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
    %% PLOT DENSITIES
    if ft == 3 && beta < 512 && (loopbeta >= 50 || change <= 0.01)
    beta = 2*beta;
    loopbeta = 0;
    change = 1;
    fprintf('Parameter beta increased to %g.\n',beta);
  end

    
end
figure('Position', [100, 100, 1200, 500]);  % 增加图像的宽度，将整张图像拉长

% 调整结构图像的大小和位置
subplot('Position', [0.05, 0.15, 0.4, 0.7]);  % 调整位置和大小 [左, 下, 宽, 高]
colormap(gray);
imagesc(1-xPhys); 
caxis([0 1]); 
axis equal; 
axis off; 
title('Optimized Structure');

% 调整目标函数历史图像的大小和位置
subplot('Position', [0.55, 0.15, 0.4, 0.7]);  % 调整位置和大小 [左, 下, 宽, 高]
plot(1:loop, c_history(1, :), '-o', 'Color', 'r', 'LineWidth', 2);  % 红色，圆形标记
hold on;
plot(1:loop, c_history(2, :), '-s', 'Color', 'g', 'LineWidth', 2);  % 绿色，方形标记
plot(1:loop, c_history(3, :), '-d', 'Color', 'b', 'LineWidth', 2);  % 蓝色，菱形标记
xlabel('Iteration');
ylabel('Objective Function Value');
title('Objective Function History');
legend('Load case 1', 'Load case 2', 'Load case 3');  % 添加图例以便区分

% 保存图像
save_filename = sprintf('MinMax_iter%d_obj%.4f.png', loop, c);
saveas(gcf, save_filename);
