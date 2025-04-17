%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE Nov, 2010 %%%%
% MMA + density + sensitivity filter (通过平凡约束实现无约束)
% discretization (nelx*nely)
% penal: the penalization power
% rmin: filter size
clear;close all;clc;
nelx=240;
nely=80;
volfrac=1e9;       % 关键修改：将体积分数设为极大值使约束失效
penal=3.0;
rmin=9.6;
ft=2;
%% MATERIAL PROPERTIES
E0 = 1;      % Young modulus of solid
Emin = 1e-9; 
nu = 0.3;    
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
F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);
fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
alldofs = [1:2*(nely+1)*(nelx+1)];
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
xval = repmat(0.1,nely,nelx); % 初始密度设为中间值
xPhys = xval;
if ft == 2
    xPhys( : ) = reshape( ( H * xval( : ) ) ./ Hs, nely, nelx );
end
loop = 0;
change = 1;
c_history = [];
%% INITIALIZE MMA OPTIMIZER
m     = 1;                % 仍然保留一个约束
n     = nelx*nely;        
xmin  = zeros(n,1);    
xmax  = ones(n,1);     
xold1 = xval(:);          
xold2 = xval(:);          
low   = [];             
upp   = [];             
a0    = 1;                
a     = zeros(m,1);       
c_MMA = 1e-3*ones(m,1);   % 关键修改：降低约束惩罚系数
d     = zeros(m,1);       
move  = 0.2;              
tic;
%% START ITERATION
while change > 0.01
    loop = loop + 1;
    %% FE-ANALYSIS
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
    K  = sparse(iK,jK,sK); K = (K+K')/2;
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);            
    c  = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));                      
    dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;                            
    dv = ones(nely,nelx);  % 保留灵敏度计算但会被大幅缩放
    %% FILTERING/MODIFICATION OF SENSITIVITIES
    if ft == 2
        dc(:) = H*(dc(:)./Hs);
        dv(:) = H*(dv(:)./Hs); % 保留密度过滤
    end
    %% METHOD OF MOVING ASYMPTOTES
    xmin = max(xval(:) - move, 0);                 
    xmax = min(xval(:) + move, 1);
    f0val = c;                                    
    df0dx = dc(:);                                
    % 关键修改：约束计算（自动产生极小值）
    fval  = sum(xPhys(:))/(volfrac*n) - 1;       % 分子最大为n, 计算结果约-1
    dfdx  = dv(:)'/(volfrac*n);                  % 梯度约为1e-9量级
    [xmma, ~, ~, ~, ~, ~, ~, ~, ~, low, upp] = ...
        mmasub(m, n, loop, xval(:), xmin, xmax, xold1, xold2,...
        f0val, df0dx, fval, dfdx, low, upp, a0, a, c_MMA, d);
    %% UPDATE VARIABLES
    xold2 = xold1;
    xold1 = xval(:);
    xval  = reshape(xmma,nely,nelx);
    %% FILTERING
    if ft == 2
        xPhys(:) = (H*xval(:))./Hs;
    else
        xPhys = xval;
    end        
    change = max(abs(xval(:)-xold1(:)));
    %% PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f ch.:%7.3f\n',loop,c,change);
    %% PLOT DENSITIES
    colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
end
time=toc;
%% FINAL VISUALIZATION
figure('Position', [100, 100, 1200, 500]);
subplot('Position', [0.05, 0.15, 0.4, 0.7]);
colormap(gray);
imagesc(1-xPhys); 
caxis([0 1]); 
axis equal; 
axis off; 
title('Optimized Structure');
subplot('Position', [0.55, 0.15, 0.4, 0.7]);
plot(1:loop, c_history, '-o', 'LineWidth', 2);
xlabel('Iteration');
ylabel('Objective Function Value');
title('Objective Function History');
saveas(gcf, sprintf('TrivialConstraint_nelx%d_time%.4f.png', nelx, time));
close all;