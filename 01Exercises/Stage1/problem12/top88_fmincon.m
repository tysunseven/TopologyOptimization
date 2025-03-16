function top88_fmincon(nelx,nely,volfrac,penal,rmin)
global c_history; % 声明全局变量
c_history = [];   % 初始化目标函数历史记录
%% MATERIAL PROPERTIES
E0 = 1;
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
%% STORE FILTER AS A SINGLE SPARSE MATRIX
HsH = spdiags(1.0./Hs,0,nelx*nely,nelx*nely)*H;
%% PACK DATA IN ONE STRUCTURE FOR USE IN SUBROUTINES
data.KE = KE; data.iK = iK; data.jK = jK; 
data.nelx = nelx; data.nely = nely; data.F = F;
data.freedofs = freedofs; data.edofMat = edofMat;
data.penal = penal; 
data.Emin = Emin; data.E0 = E0; data.volfrac = volfrac;
data.HsH  = HsH;
%% VOLUME CONSTRAINT
A = (data.HsH'*repmat(1/(nelx*nely),nelx*nely,1))';
%% LOWER AND UPPER DESIGN BOUNDS
xmin = zeros(nelx*nely,1);
xmax = ones (nelx*nely,1);
%% INITIAL DESIGN
x0 = volfrac(ones(nelx*nely,1));
%% OPTIMIZATION OPTIONS (help fmincon)
options = optimset('Algorithm', 'interior-point', ...
                   'GradObj','on', ... % 'On' indicates that the user needs to provide the gradient themselves.
                   'DerivativeCheck', 'off', ...
                   'PlotFcn', {@(x,ov,s)PlotX(x,ov,s,data)}, ...                   
                   'Hessian',{'lbfgs',25}, ...
                   'MaxFunEvals', 10000, ...  % Corresponds to MaxFunctionEvaluations
                   'MaxIter', 1000, ...  % Corresponds to MaxIterations
                   'OutputFcn', @myOutputFcn);  % Set constraint tolerance);
%'Display','iter', ...
%'TolFun', 1e6, ...  % Corresponds to FunctionTolerance
 %                  'TolCon', 1e6, ...
%% CALL THE OPTIMIZATION SOFTWARE
tic;
[xFinal,fval,flag,output] = fmincon( ...
    ... % objective function/gradient:
    @(x) fdf(x,data), ...
    ... % initial guess:
    x0,...
    ... % linear inequality constraints: A*x <= volfrac
    A, volfrac, ...
    ... % linear equality constraints: none
    [],[], ...
    ... % lower/upper bounds
    xmin, xmax, ...
    ... % non-linear constraints: none
    [], ...
    ... % finally, optimization options
    options);

xFinal = data.HsH * xFinal;  % 计算最终设计变量的物理密度

time = toc;
% 使用 length(c_history) 获取迭代次数
figure('Position', [100, 100, 1200, 500]);  % 增加图像的宽度，将整张图像拉长

% 调整结构图像的大小和位置
subplot('Position', [0.05, 0.15, 0.4, 0.7]);  % 调整位置和大小 [左, 下, 宽, 高]
colormap(gray);
imagesc(1-reshape(xFinal,data.nely,data.nelx)); 
caxis([0 1]); 
axis equal; 
axis off; 
title('Optimized Structure');

% 调整目标函数历史图像的大小和位置
subplot('Position', [0.55, 0.15, 0.4, 0.7]);  % 调整位置和大小 [左, 下, 宽, 高]
plot(1:length(c_history), c_history, '-o', 'LineWidth', 2);  % 使用 -o 添加数据点标记
xlabel('Iteration');
ylabel('Objective Function Value');
title('Objective Function History');

% 保存图像
save_filename = sprintf('IP_nelx%d_time%.4f_iter%d_obj%.4f.png', nelx, time, length(c_history), fval);
saveas(gcf, save_filename);
close all;


%------------------------------------------------------------------
%% Objective function
function [c,dc] = fdf(x,data)
if(nargout>1)
  [c,dc]=FE(x,data);
else
  [c]=FE(x,data);
end
%------------------------------------------------------------------
%% Plot function
function [stop]=PlotX(x,~,~,data)
%% Filtering
xPhys = data.HsH * x; 
%% PLOT DENSITIES
colormap(gray); imagesc(1-reshape(xPhys,data.nely,data.nelx)); 
clim([0 1]); axis ij; axis equal; axis off; drawnow;
stop=false;
%------------------------------------------------------------------
%% perform FE-analysis
function [c,dc]=FE(x,d)
% store一些变量 between the calls
global c_history;  % 声明全局变量 c_history
persistent U x_old xPhys L s ce
if length(x_old) ~= length(x)
  % pre-allocate memory
  x_old = repmat(NaN,size(x));
  U     = zeros(2*(d.nely+1)*(d.nelx+1),1); 
  % A     = zeros(2*(d.nely+1)*(d.nelx+1),1); 
end
%
if any(x_old~=x)
  % need to re-assemble and re-factorize
  x_old = x;
  %% Filtering
  xPhys = d.HsH * x;
  %% FE-ANALYSIS
  sK = reshape(d.KE(:)*(d.Emin+xPhys(:)'.^d.penal*(d.E0-d.Emin)),...
               64*d.nelx*d.nely,1);
  K = sparse(d.iK,d.jK,sK); K = (K+K')/2;
  % Cholesky factorization
  [L,~,s]=chol(K(d.freedofs,d.freedofs),'lower','vector');
  % Forward/backward substitution
  U(d.freedofs(s))=L'\(L\d.F(d.freedofs(s)));
  %
  ce =  sum((U(d.edofMat)*d.KE).*U(d.edofMat),2);
end
%
% compute outputs
%
%FIXIT: compute objective function (compliance)
c = sum((d.Emin+xPhys(:).^d.penal).*(d.E0-d.Emin).*ce);
if nargout > 1
    %FIXIT: compute sensitivities of compliance
    dc = -d.penal*(d.E0-d.Emin)*xPhys.^(d.penal-1).*ce;
    %% MODIFICATION OF SENSITIVITIES
    dc = d.HsH' * dc;
end


function stop = myOutputFcn(x, optimValues, state)
    persistent x_old1  % 保存上一轮的设计变量
    global c_history;  % 声明全局变量 c_history
    
    stop = false;  % 默认情况下不终止优化

    switch state
        case 'init'
            disp('Starting the optimization.');
            x_old1 = x;  % 初始化 x_old 为第一轮设计变量
        case 'iter'
            % 计算当前设计变量与上一轮设计变量的最大变化量
            max_change = max(abs(x - x_old1));
            
            % 计算目标函数值 (取自 optimValues.fval)
            c = optimValues.fval;  % 当前的目标函数值

            % 将当前的目标函数值添加到 c_history 中
            c_history = [c_history; c];
            
            % 打印与目标类似的输出格式
            fprintf(' It.:%5i Obj.:%11.4f  ch.:%7.3f\n', ...
                optimValues.iteration, c, max_change);

            if optimValues.iteration > 10 && max_change < 0.01
                disp('Terminating optimization: max_change is smaller than tolX.');
                stop = true;  % 手动终止优化
            end
            
            % 更新上一轮的设计变量为当前轮的设计变量
            x_old1 = x;
        case 'done'
            disp('Optimization finished.');
        otherwise
            % 处理其他状态
    end