%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE Nov, 2010 %%%%
% MMA + density + sensitivity filter
% discretization (nelx*nely)
% penal: the penalization power
% rmin: filter size
clear;close all;clc;
nelx=240;
nely=80;
volfrac=0.5;
penal=3.0;  %  typically 3
rmin=9.6;
ft=2;       %  ft=1 is sensitivity filter, ft=2 is density filter
%% MATERIAL PROPERTIES
E0 = 1;      % Young modulus of solid
Emin = 1e-9; % Young modulus of " void "
nu = 0.3;    % Poisson ratio
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);                          % nodes of squre mesh for entire domain
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);                   % first dof in every elemnts
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);% each row is the 8 dofs for every elemnts
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
xval = repmat(volfrac,nely,nelx);
xPhys = xval;
if ft == 1                                                                 
    xPhys = xval;
elseif ft == 2
    xPhys( : ) = reshape( ( H * xval( : ) ) ./ Hs, nely, nelx );
end
loop = 0;
change = 1;
c_history = [];
%% INITIALIZE MMA OPTIMIZER
% DESIGN UPDATE BY THE MMA METHOD
m     = 1;                % The number of general constraints. Volumn constrain
n     = nelx*nely;        % The number of design variables x_j.
xmin  = zeros( n, 1 );    % Column vector with the lower bounds for the variables x_j.
xmax  = ones( n, 1 );     % Column vector with the upper bounds for the variables x_j.
xold1 = xval(:);          % xval, one iteration ago (provided that iter>1).
xold2 = xval(:);          % xval, two iterations ago (provided that iter>2).
low   = [  ];             % Column vector with the lower asymptotes from the previous iteration (provided that iter>1).
upp   = [  ];             % Column vector with the upper asymptotes from the previous iteration (provided that iter>1).
a0    = 1;                % The constants a_0 in the term a_0*z.
a     = zeros(m,1);       % Column vector with the constants a_i in the terms a_i*z.
c_MMA = 1000*ones(m,1);   % Column vector with the constants c_i in the terms c_i*y_i.
d     = zeros(m,1);       % Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
move  = 0.2;              % Max allowed change in design variables
tic;
%% START ITERATION
while change > 0.01
    loop = loop + 1;
    %% FE-ANALYSIS
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
    K  = sparse(iK,jK,sK); K = (K+K')/2;
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);            % ce=U'*K*U
    c  = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));                      % eq 2 in top88
    dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;                            % eq 5 in top88
    dv = ones(nely,nelx);                                                  % eq 6 in top88
    c_history = [c_history c];  % Store objective function value
    %% FILTERING/MODIFICATION OF SENSITIVITIES
    %  ft=1 is sensitivity filter, ft=2 is density filter
    if ft == 1
        dc(:) = H*(xval(:).*dc(:))./Hs./max(1e-3,xval(:));
    elseif ft == 2
        dc(:) = H*(dc(:)./Hs);
        dv(:) = H*(dv(:)./Hs);
    end
%% METHOD OF MOVING ASYMPTOTES
    xmin = max( xval(:) - move, 0 );                 
    xmax = min( xval(:) + move, 1 );
    f0val = c;                                    % scalar
    df0dx = dc(:);                                % column 
    fval  = sum( xPhys(:))/(volfrac*n) - 1;       % scalar volumn constrain
    dfdx  = dv(:)'/ (volfrac*n);                  % row
    [ xmma, ~, ~, ~, ~, ~, ~, ~, ~, low, upp ] = ...
        mmasub( m, n, loop, xval(:), xmin, xmax, xold1, xold2, f0val,df0dx,fval,dfdx,low,upp,a0,a,c_MMA,d);
    % Update MMA Variables
    xold2    = xold1(:);
    xold1    = xval(:);
    xval     = reshape(xmma,nely,nelx);
    % FILTERING
    if ft == 1         % sensitivity filter
        xPhys = xval;
    elseif ft == 2     % density filter
        xPhys(:) = (H *xval(:))./ Hs;
    end        
    change = max(abs(xval(:)-xold1(:)));
      
    %% PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, ...
        mean(xPhys(:)),change);
    %% PLOT DENSITIES
    colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
   
end
time=toc;
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
plot(1:loop, c_history, '-o', 'LineWidth', 2);  % 使用 -o 添加数据点标记
xlabel('Iteration');
ylabel('Objective Function Value');
title('Objective Function History');

% 保存图像
save_filename = sprintf('MMA_nelx%d_time%.4f_iter%d_obj%.4f.png', nelx, time, loop, c);
saveas(gcf, save_filename);
close all;