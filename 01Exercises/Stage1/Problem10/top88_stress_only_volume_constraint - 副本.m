%% Set Parameters

nelx=40;
nely=20;
volfrac=0.3;
penal=3.0;
rmin=nely*0.06;
ft=3;
beta = 1; % 初始Heaviside参数
pnorm = 12;
normalized_coefficient = 1;
vmstress_pnorm = 1;


%% MATERIAL PROPERTIES

E0 = 100;
Emin = 1e-9*E0;
nu = 0.3;
In_spring=0.1;
Out_spring=0.1;

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

din = 1; dout = 2*(nely+1)*nelx+1;
fixeddofs =union([2*(nely+1)-1,2*(nely+1)],2:2*(nely+1):2*(nely+1)*nelx+2);%fix x of top left and y of bottom edge

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
vm_stress = zeros(nely,nelx);
relaxed_stress = zeros(nely, nelx);
xPhys = xval;
xTilde = xval;
if ft == 1
    xPhys = xval;
elseif ft == 2
    xPhys( : ) = reshape( ( H * xval( : ) ) ./ Hs, nely, nelx );
elseif ft == 3
    xPhys = 1 - exp(-beta * xTilde) + xTilde * exp(-beta);
end
loop = 0;
loopbeta = 0;
c_history = [];
change = 1;

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
fval = zeros(m,1);
dfdx = zeros(m,n);


figure;

%% START ITERATION

while change > 0.003
    loop = loop + 1;
    loopbeta = loopbeta+1;

    if loop > 1
        normalized_coefficient = 0.5 * max(relaxed_stress(:))/vmstress_pnorm+0.5*normalized_coefficient;
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

    

    vm_stress = double(reshape(sqrt(max(sum((U(edof) * Q) .* U(edof), 2), 0)), nely, nelx));
    relaxed_stress = xPhys.^2 .* vm_stress;
    vmstress_p = relaxed_stress.^pnorm;
    sum_vmstress_p = sum(vmstress_p(:));
    vmstress_pnorm = sum_vmstress_p^(1/pnorm);
    normalized_vmstress_pnorm = normalized_coefficient * vmstress_pnorm;
    % normalized_vmstress_pnorm = normalized_coefficient * sum(vmstress(:).^pnorm)^(1/pnorm);

    dfsdu = normalized_coefficient * sum_vmstress_p^(1/pnorm-1);

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
    f0val = c;
    df0dx = dc(:);
    fval(1)  = sum(xPhys(:))/(volfrac*n) - 1;
    dfdx(1,:)  = dv(:)'/ (volfrac*n);
    % fval(2)  = sum(xPhys(:))/(volfrac*n) - 1;
    % dfdx(2,:)  = dv(:)'/ (volfrac*n);
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

    %% PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, ...
        mean(xPhys(:)),change);

    subplot(1, 2, 1);  % 创建 1 行 2 列的子图，并选择第 1 个位置
    colormap(gray);
    imagesc(1-xPhys);  % 绘制矩阵图像
    title('1 - xPhys');  % 添加标题
    axis equal;  % 设置坐标轴比例相同
    axis off;

    subplot(1, 2, 2);  % 选择第 2 个子图位置
    colormap(gray);
    imagesc(-relaxed_stress);  % 绘制应力矩阵图像
    title('vmstress');  % 添加标题
    axis equal;  % 设置坐标轴比例相同
    axis off;



    drawnow;
end
% save_filename = sprintf('Without_stress_constraints_iter%d_obj%.4f_maxstress%.4f.png', loop, c, max(vmstress(:)));
% saveas(gcf, save_filename);





