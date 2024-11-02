%% Set Parameters
nelx=50;
nely=25;
nele = nely * nelx;
ndof = (nely + 1) * (nelx + 1);
change_threshold = 0.001;
volfrac=0.3;
penal=3.0;
rmin=nely*0.06;
ft=3;
m  = 1;                % The number of general constraints. Volumn constrain

beta = 1;
p = 12;
q = 2;
sigma = 1.15;


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
Q_vector = Q(:);

normalized_coefficient = 1;


%% PREPARE FINITE ELEMENT ANALYSIS

A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);

% 第i行存储了第i个单元（对单元的排序是先从上到下遍历第一列，再第二列...）的8个自由度在总的自由度排序中的索引
% 对某个单元的8个自由度排序是从左下节点的x分量出发，逆时针转
% edofMat的元素的取值范围是1到2*ndof
% 如果有一个关于节点的量，比如U，它是一个2*nodf维的列向量，那U(edofMat)就可以提取出得到每个单元的这个量的信息
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);

% 考虑一个8*8的矩阵,将单元按照我们惯用的顺序排序
% 从上到下数单元的横坐标，就是8个12345678，得到一个64行的列向量，不妨记作I
% 从上到下数单元的纵坐标，就是8个1，8个2，...，8个8，得到一个64行的列向量，不妨记作J
% 在实践中我们会将这两个列向量与edofMat的第i行复合，得到两个64行的列向量edofMat(i,I)和edofMat(i,J)
% 在将某个局部矩阵拼到某个整体矩阵的时候，同时遍历edofMat(i,I)、edofMat(i,J)和局部矩阵展开成的列向量，就知道该把局部矩阵的元素拼到哪里
% 把edofMat(i,I)和edofMat(i,J)按照i的顺序拼接成更长的列向量，第i个64行就记录了该怎么从局部拼到整体。重新将这两个64*nelx*nely行的列向量记为I和J
% 本题中有第i个局部矩阵v_iQ需要被拼到整体矩阵，经过分析，拼法就是上面的拼法
I = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
J = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);


%% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)

U = zeros(2*(nely+1)*(nelx+1),1);
Lambda_objective = zeros(2*(nely+1)*(nelx+1),1);
Lambda_stress = zeros(2*(nely+1)*(nelx+1),1);
L_out= zeros(2*(nely+1)*(nelx+1),1);

din = 1; dout = 2*(nely+1)*nelx+1;
fixeddofs =union([2*(nely+1)-1,2*(nely+1)],2:2*(nely+1):2*(nely+1)*nelx+2);%fix x of top left and y of bottom edge

F = sparse(din,1,1,2*(nely+1)*(nelx+1),1);
alldofs = 1:2*(nely+1)*(nelx+1);
freedofs = setdiff(alldofs,fixeddofs);

fixedeles = union(1,nely*(nelx-1)+1);
alleles= 1:nely*nelx;
freeeles = setdiff(alleles,fixedeles);

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
dfsdxPhysPart1 = zeros(1,n);
dfsdxPhysPart2 = zeros(1,n);
dfsdxPhys = zeros(1,n);

figure;

%% START ITERATION

while change > change_threshold && m==1
    loop = loop + 1;
    loopbeta = loopbeta+1;


    % 算K
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
    K = sparse(I,J,sK); K = (K+K')/2; K(din,din) = K(din,din) + In_spring;  K(dout,dout) = K(dout,dout) + Out_spring;
    % 算U
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);
    % 算目标函数
    c = U(dout);
    c_history = [c_history c];
    % 伴随方法算目标函数的灵敏度
    L_out(dout)=1;
    Lambda_objective(freedofs) = -K(freedofs,freedofs)\L_out(freedofs);
    ce = reshape(sum((Lambda_objective(edofMat)*KE).*U(edofMat),2),nely,nelx);
    dc = penal*(E0-Emin)*xPhys.^(penal-1).*ce;
    dv = ones(nely,nelx);

    Ue = U(edofMat);
    vm_stress = double(reshape(sqrt(max(sum((Ue * Q) .* Ue, 2), 0)), nely, nelx));
    relaxed_stress = xPhys.^q .* vm_stress;
    sum_relaxed_stress_p = sum(relaxed_stress(freeeles).^p);
    fs = normalized_coefficient * sum_relaxed_stress_p^(1/p);

    part1 = sum_relaxed_stress_p^(1/p-1) * q * reshape(xPhys.^(q*p-1) .* vm_stress.^p, nele, 1);
    part2V = reshape(xPhys.^(q*p) .* vm_stress.^(p-2), nele, 1);
    part2V(fixedeles) = 0;
    spart2 = reshape(Q_vector*part2V',nele*64,1);
    S = sparse(I, J, spart2, 2*ndof, 2*ndof);
    SU = S*U;

    Lambda_stress(freedofs) = -K(freedofs,freedofs)\ SU(freedofs);
    part2 = sum_relaxed_stress_p^(1/p-1) * penal*(E0-Emin)*xPhys(:).^(penal-1).*sum((Lambda_stress(edofMat)*KE).*U(edofMat),2);
    dfsdxPhys = normalized_coefficient*(part1 + part2);
    dfsdxPhys(fixedeles) = 0;

    
    normalized_coefficient = 0.5 * max(relaxed_stress(:))/sum_relaxed_stress_p^(1/p)+0.5*normalized_coefficient;
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
        dfsdxPhys = H*(dfsdxPhys.*dx(:)./Hs);
    end

    %% METHOD OF MOVING ASYMPTOTES

    xmin = max( xval(:) - move, 0 );
    xmax = min( xval(:) + move, 1 );
    f0val = c;
    df0dx = dc(:);
    fval(1)  = sum(xPhys(:))/(volfrac*n) - 1;
    dfdx(1,:)  = dv(:)'/ (volfrac*n);

    if m == 2
    % fval(2)  = sum(xPhys(:))/(volfrac*n) - 1;
    fval(2) = fs - sigma;
    % dfdx(2,:)  = dv(:)'/ (volfrac*n);
    dfdx(2,:)  = dfsdxPhys';
    end;

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
    if beta<4
        if ft == 3 && beta < 512 && (loopbeta >= 51 || change <= 2*change_threshold)
            beta = 2*beta;
            loopbeta = 0;
            change = 1;
            fprintf('Parameter beta increased to %g.\n',beta);
        end
    end
    change = max(abs(xval(:)-xold1(:)));

    %% PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f stress.:%7.3f beta.:%5i\n',loop,c, ...
        mean(xPhys(:)),change,max(relaxed_stress(freeeles)),beta);

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
    if loop>75
        m=2;
        loop=0;
        loopbeta=0;
    end
end




if change<change_threshold
    return
end

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
dfsdxPhysPart1 = zeros(1,n);
dfsdxPhysPart2 = zeros(1,n);
dfsdxPhys = zeros(1,n);

figure;

%% START ITERATION

while change > change_threshold
    loop = loop + 1;
    loopbeta = loopbeta+1;


    % 算K
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
    K = sparse(I,J,sK); K = (K+K')/2; K(din,din) = K(din,din) + In_spring;  K(dout,dout) = K(dout,dout) + Out_spring;
    % 算U
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);
    % 算目标函数
    c = U(dout);
    c_history = [c_history c];
    % 伴随方法算目标函数的灵敏度
    L_out(dout)=1;
    Lambda_objective(freedofs) = -K(freedofs,freedofs)\L_out(freedofs);
    ce = reshape(sum((Lambda_objective(edofMat)*KE).*U(edofMat),2),nely,nelx);
    dc = penal*(E0-Emin)*xPhys.^(penal-1).*ce;
    dv = ones(nely,nelx);

    Ue = U(edofMat);
    vm_stress = double(reshape(sqrt(max(sum((Ue * Q) .* Ue, 2), 0)), nely, nelx));
    relaxed_stress = xPhys.^q .* vm_stress;
    sum_relaxed_stress_p = sum(relaxed_stress(freeeles).^p);
    fs = normalized_coefficient * sum_relaxed_stress_p^(1/p);

    part1 = sum_relaxed_stress_p^(1/p-1) * q * reshape(xPhys.^(q*p-1) .* vm_stress.^p, nele, 1);
    part2V = reshape(xPhys.^(q*p) .* vm_stress.^(p-2), nele, 1);
    part2V(fixedeles) = 0;
    spart2 = reshape(Q_vector*part2V',nele*64,1);
    S = sparse(I, J, spart2, 2*ndof, 2*ndof);
    SU = S*U;

    Lambda_stress(freedofs) = -K(freedofs,freedofs)\ SU(freedofs);
    part2 = sum_relaxed_stress_p^(1/p-1) * penal*(E0-Emin)*xPhys(:).^(penal-1).*sum((Lambda_stress(edofMat)*KE).*U(edofMat),2);
    dfsdxPhys = normalized_coefficient*(part1 + part2);
    dfsdxPhys(fixedeles) = 0;

    
    normalized_coefficient = 0.5 * max(relaxed_stress(:))/sum_relaxed_stress_p^(1/p)+0.5*normalized_coefficient;
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
        dfsdxPhys = H*(dfsdxPhys.*dx(:)./Hs);
    end

    %% METHOD OF MOVING ASYMPTOTES

    xmin = max( xval(:) - move, 0 );
    xmax = min( xval(:) + move, 1 );
    f0val = c;
    df0dx = dc(:);
    fval(1)  = sum(xPhys(:))/(volfrac*n) - 1;
    dfdx(1,:)  = dv(:)'/ (volfrac*n);


    % fval(2)  = sum(xPhys(:))/(volfrac*n) - 1;
    fval(2) = fs - sigma;
    % dfdx(2,:)  = dv(:)'/ (volfrac*n);
    dfdx(2,:)  = dfsdxPhys';


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

    if ft == 3 && beta < 512 && (loopbeta >= 50 || change <= 2*change_threshold)
        beta = 2*beta;
        loopbeta = 0;
        change = 1;
        fprintf('Parameter beta increased to %g.\n',beta);
    end

    change = max(abs(xval(:)-xold1(:)));

    %% PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f stress.:%7.3f beta.:%5i\n',loop+75,c, ...
        mean(xPhys(:)),change,max(relaxed_stress(freeeles)),beta);

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

save_filename = sprintf('With_stress_constraints_iter%d_obj%.4f_stress%.4f.png', loop+75, c, max(relaxed_stress(freeeles)));
saveas(gcf, save_filename);