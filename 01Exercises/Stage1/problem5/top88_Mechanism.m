%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE Nov, 2010 %%%%
%discretization (nelx*nely)
clear;close all;clc;
nelx=40;
nely=40;
volfrac=0.3;
penal=3.0;      % the penalization power, typically 3
rmin=2;         % filter size
ft=3;
%% MATERIAL PROPERTIES
E0 = 100;       % Young modulus of solid
Emin = 1e-9*E0; % Young modulus of " void "
nu = 0.3;       % Poisson ratio
In_spring=0.1;
Out_spring=0.1;
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
%% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
U = zeros(2*(nely+1)*(nelx+1),1);
Lambda = zeros(2*(nely+1)*(nelx+1),1); % compliant add to this
L= zeros(2*(nely+1)*(nelx+1),1);       % compliant add to this

% BC with symmetry
% din = (nely+1)*2-1; dout = 2*(nely+1)*(nelx+1)-1;
% fixeddofs =union([1,2],2*(nely+1):2*(nely+1):2*(nely+1)*(nelx+1));%fix x of top left and y of bottom edge

% BC without symmetry
din=2*ceil((nely+1)/2)-1;dout = 2*(nelx)*(nely+1)+2*ceil((nely+1)/2)-1;
fixeddofs = union([1:2], [2*(nely+1)-1,2*(nely+1)]);

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
change = 1;
loopbeta = 0;
c_history = [];
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
%% START ITERATION
while change > 0.001
    loopbeta = loopbeta+1;
    loop = loop + 1;
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
    %% FILTERING/MODIFICATION OF SENSITIVITIES
    %  ft=1 is sensitivity filter, ft=2 is density filter
    if ft == 1
        dc(:) = H*(xval(:).*dc(:))./Hs./max(1e-3,xval(:));
    elseif ft == 2
        dc(:) = H*(dc(:)./Hs);
        dv(:) = H*(dv(:)./Hs);
    elseif ft == 3 % Heaviside投影的敏感度滤波
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
    if ft == 1         %sensitivity filter
        xPhys = xval;
    elseif ft == 2     %density filter
        xPhys(:) = (H *xval(:))./ Hs;
        elseif ft == 3
    xTilde(:) = (H * xval(:))./ Hs;
    xPhys = 1 - exp(-beta * xTilde) + xTilde * exp(-beta); % Heaviside正则化
    end
    if ft == 3 && beta < 512 && (loopbeta >= 50 || change <= 0.01)
    beta = 2*beta;
    loopbeta = 0;
    change = 1;
    fprintf('Parameter beta increased to %g.\n',beta);
  end
    change = max(abs(xval(:)-xold1(:)));

    %% PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, ...
        mean(xPhys(:)),change);
    colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
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
plot(1:loop, c_history, '-o', 'LineWidth', 2);  % 使用 -o 添加数据点标记
xlabel('Iteration');
ylabel('Objective Function Value');
title('Objective Function History');

% 保存图像
save_filename = sprintf('Mechanism_ft%d_iter%d_obj%.4f.png', ft, loop, c);
saveas(gcf, save_filename);
close all;
        