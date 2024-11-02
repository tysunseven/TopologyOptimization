%% Set Parameters

nelx=40;
nely=40;
volfrac=0.3;
penal=3.0;
rmin=2;

ft=3;
gamma=1;

%% MATERIAL PROPERTIES

E0 = 100;
Emin = 1e-9*E0;
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

%% DEFINE LOADS AND SUPPORTS (COMPLIANT INVERTER MECHANISM)

U = zeros(2*(nely+1)*(nelx+1),3);
Lambda = zeros(2*(nely+1)*(nelx+1),3);

din = 2*ceil((nely+1)/2)-1;
dout = 2*(nelx)*(nely+1)+2*ceil((nely+1)/2)-1;
In_spring = 0.1;
Out_spring = 0.1;

F = sparse(din,1,1,2*(nely+1)*(nelx+1),1);
L = zeros(2*(nely+1)*(nelx+1),1);
L(dout)=1;


alldofs = 1:2*(nely+1)*(nelx+1);
fixeddofs = union([1:2], [2*(nely+1)-1,2*(nely+1)]);
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

beta = 1;
if ft == 1 || ft == 2
    xPhys(:,:,1) = xval;
    xPhys(:,:,2) = xval;
    xPhys(:,:,3) = xval;
elseif ft == 3
    xTilde = xval;
    xPhys = zeros(nely, nelx, 3);
    xPhys(:,:,1) = 1-exp(-beta*xTilde)+xTilde*exp(-beta);
    xPhys(:,:,2) = 1-exp(-gamma*xPhys(:,:,1))+xPhys(:,:,1)*exp(-gamma);
    xPhys(:,:,3) = exp(-gamma*(1-xPhys(:,:,1)))-(1-xPhys(:,:,1))*exp(-gamma);
end

loop = 0;
loopbeta = 0;

change = 1;

c_history = [];

%% INITIALIZE MMA OPTIMIZER

p=3;                      % ideal+dilate+erode
q=1;                      % volume constrain
m = p+q;                % The number of general constraints.
n = nelx*nely;            % The number of design variables.
xmin  = zeros(n,1);       % Column vector with the lower bounds for the variables x_j.
xmax  = ones(n,1);        % Column vector with the upper bounds for the variables x_j.
xold1 = xval(:);          % xval, one iteration ago (provided that iter>1).
xold2 = xval(:);          % xval, two iterations ago (provided that iter>2).
low   = [  ];             % Column vector with the lower asymptotes from the previous iteration (provided that iter>1).
upp   = [  ];             % Column vector with the upper asymptotes from the previous iteration (provided that iter>1).
a0    = 1;                % The constants a_0 in the term a_0*z.
a     = 1*ones(m,1);      % a_i = 1, i=1,...,2p
a(end)=0;                 % a_{2p+i} = 0, i=q
c_MMA = 1000*ones(m,1);   % Column vector with the constants c_i in the terms c_i*y_i.
d     = zeros(m,1);       % Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
move  = 0.1;              % Max allowed change in design variables
fval = zeros(m,1);
dfdx = zeros(m,n);

% 创建图形窗口
figure('Position', [100, 100, 1200, 500]);  % 增加图像的宽度
subplot(1, 3, 1);  % 1行3列的第1个子图
colormap(gray);
title('Optimized Structure 1');
axis equal;
axis off;

subplot(1, 3, 2);  % 1行3列的第2个子图
title('Optimized Structure 2');
axis equal;
axis off;

subplot(1, 3, 3);  % 1行3列的第3个子图
title('Optimized Structure 3');
axis equal;
axis off;

%% START ITERATION
while change > 0.002

    loop = loop + 1;
    loopbeta = loopbeta+1;

    %% FE-ANALYSIS


    xPhys1=xPhys(:,:,1);
    xPhys2=xPhys(:,:,2);
    xPhys3=xPhys(:,:,3);

    sk=zeros(64*nelx*nely,3);
    sK(:,1) = reshape(KE(:)*(Emin+xPhys1(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
    sK(:,2) = reshape(KE(:)*(Emin+xPhys2(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
    sK(:,3) = reshape(KE(:)*(Emin+xPhys3(:)'.^penal*(E0-Emin)),64*nelx*nely,1);

    K = cell(1, 3);
    for i = 1:3
        K{i} = sparse(iK, jK, sK(:, i));
        K{i} = (K{i} + K{i}') / 2;
        K{i}(din, din) = K{i}(din, din) + In_spring;
        K{i}(dout, dout) = K{i}(dout, dout) + Out_spring;
        Lambda(freedofs,i) = -K{i}(freedofs,freedofs)\L(freedofs);
    end

    for i = 1:3
        U(freedofs, i) = K{i}(freedofs, freedofs) \ F(freedofs);
    end

    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS

    U1=U(:,1);U2=U(:,2);U3=U(:,3);
    Lambda1=Lambda(:,1);Lambda2=Lambda(:,2);Lambda3=Lambda(:,3);
    
    c = zeros(3,1);
    dc = zeros(nely,nelx,3);

    ce1 = reshape(sum((Lambda1(edofMat)*KE).*U1(edofMat),2),nely,nelx);
    ce2 = reshape(sum((Lambda2(edofMat)*KE).*U2(edofMat),2),nely,nelx);
    ce3 = reshape(sum((Lambda1(edofMat)*KE).*U3(edofMat),2),nely,nelx);
    c(1) = U1(dout);
    c(2) = U2(dout);
    c(3) = U3(dout);
    dc(:,:,1) =  penal*(E0-Emin)*xPhys(:,:,1).^(penal-1).*ce1;
    dc(:,:,2) =  penal*(E0-Emin)*xPhys(:,:,2).^(penal-1).*ce2;
    dc(:,:,3) =  penal*(E0-Emin)*xPhys(:,:,3).^(penal-1).*ce3;
    dc1=dc(:,:,1);
    dv = ones(nely,nelx);
    
    c_history = [c_history c];  % Store objective function value
    z=min(c);

    %% FILTERING/MODIFICATION OF SENSITIVITIES
    if ft == 1
        dc(:) = H*(xval(:).*dc(:))./Hs./max(1e-3,xval(:));
    elseif ft == 2
        for i=1:3
            dc_temp = dc(:,:,i);
            dc_temp(:) = H*(dc_temp(:)./Hs);
            dc(:,:,i) = dc_temp;
        end
        dv(:) = H*(dv(:)./Hs);
    elseif ft == 3
        dx = beta * exp(-beta * xTilde) + exp(-beta);
        dc1(:) = H*(dc1(:).*dx(:)./Hs);
        dc(:,:,1)=dc1;
        dv(:) = H*(dv(:).*dx(:)./Hs);
    end
    %% METHOD OF MOVING ASYMPTOTES
    xmin = max( xval(:) - move, 0 );
    xmax = min( xval(:) + move, 1 );
    %% Modifed for Min Max formulation
    f0val = 0;
    df0dx = zeros(n, 1);

    fval(1:p) = c(1:p);
    %fval(1:p) = c(1);
    %fval(p+1:2*p) = -c(1);

    fval(m) = sum(xPhys1(:))/(volfrac*n) - 1;
    for i=1:3
        dfdx(i,:) = reshape(dc(:,:,i),1,n);
        %dfdx(i,:) = reshape(dc(:,:,1),1,n);
        %dfdx(p+i,:) = reshape(-dc(:,:,1),1,n);
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
        xPhys1(:,:,1) = reshape((H *xval(:))./ Hs,nely,nelx);
    elseif ft == 3
        xTilde(:) = (H * xval(:))./ Hs;
        xPhys(:,:,1) = 1-exp(-beta*xTilde)+xTilde*exp(-beta);
        xPhys(:,:,2) = 1-exp(-gamma*xPhys(:,:,1))+xPhys(:,:,1)*exp(-gamma);
        xPhys(:,:,3) = exp(-gamma*(1-xPhys(:,:,1)))-(1-xPhys(:,:,1))*exp(-gamma);
    end
    change = max(abs(xval(:)-xold1(:)));

    %% PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,z, ...
        mean(xval(:)),change);

    % colormap(gray); imagesc(1-xPhys1); caxis([0 1]); axis equal; axis off; drawnow;
    % 创建新的图形窗口
% 更新子图内容
    subplot(1, 3, 1);
    imagesc(1 - xPhys(:,:,1));
    caxis([0 1]);
    axis equal;
    axis off;

    subplot(1, 3, 2);
    imagesc(1 - xPhys(:,:,2));
    caxis([0 1]);
    axis equal;
    axis off;

    subplot(1, 3, 3);
    imagesc(1 - xPhys(:,:,3));
    caxis([0 1]);
    axis equal;
    axis off;

    % 更新总标题
    sgtitle('Optimized Structures');  
    drawnow;  % 更新图形窗口

    %% PLOT DENSITIES
    if ft == 3 && beta < 512 && loopbeta >= 50
        beta = 2*beta;
        loopbeta = 0;
        change = 1;
        fprintf('Parameter beta increased to %g.\n',beta);
    end

end    % while change > 0.002

%% PLOT THE FINAL IMAGE

figure('Position', [100, 100, 1200, 500]);  % 增加图像的宽度，将整张图像拉长

% 调整结构图像的大小和位置
subplot('Position', [0.05, 0.15, 0.4, 0.7]);  % 调整位置和大小 [左, 下, 宽, 高]
colormap(gray);
imagesc(1-xPhys1);
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
