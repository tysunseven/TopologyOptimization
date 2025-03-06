% 参数设置
nelx = 20;  % x方向网格数
nely = 20;  % y方向网格数
nelem = nelx * nely; % 单元总数
nnode = (nelx+1)*(nely+1); % 节点总数

% p将是一个nnode维的列向量，u将是一个2nnode维的列向量
% Ka和Ma将是nnode×nnode维的方阵,Ks和Ms将是一个2nnode×2nnode维的方阵
% Sp将是一个2nnode行nnode列的矩阵,Su将是一个nnode行2nnode列的矩阵

fprintf('Ka和Ma应是%d阶的方阵\n', nnode);
fprintf('Ks和Ms应是%d阶的方阵\n', 2*nnode);

E = 0.38; % 固体的杨氏模量
Ev = E/10^6; % 搭配固体的虚拟的杨氏模量
nu = 0.35; % 固体的泊松比
rhos = 1190; % 固体的密度
rhosv = rhos/10^7; % 搭配固体的虚拟的密度

kappa = 141834.999; % 气体的体积模量，待确定
kappav = kappa*10^7; % 搭配气体的虚拟的体积模量，待确定
rhoa = 1.225; % 气体的密度，待确定
rhoav = rhoa*10^6; % 搭配气体的虚拟的密度，待确定
% volfrac = 1; % 体积分数，也决定了初始解
q = 1; % RAMP插值的参数

x = readmatrix('micro2.xlsx');
imagesc(1-x);
colormap(gray);
colorbar;
axis equal;
axis tight;

% x = repmat(volfrac,nely,nelx);
xPhys = x; % 暂时直接令成 x，后面再补充过滤和Heaviside投影

%% PREPARE FINITE ELEMENT ANALYSIS

nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx); % 给(nely+1)*(nelx+1)个节点编号
edofVeca = reshape(nodenrs(1:end-1,1:end-1)+nely+1,nelx*nely,1); % 记录每个单元的右上节点的编号
edofMata = repmat(edofVeca,1,4)+repmat([0 -nely-1 -nely 1],nelx*nely,1); % 记录每个单元的四个声压自由度的编号
edofVecs = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1); % 记录每个单元的左下节点的编号
edofMats = repmat(edofVecs,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1); % 记录每个单元的八个位移自由度的编号

% 用来拼Ka的单个矩阵
Kae = 1/6*[4 -1 -2 -1;
           -1 4 -1 -2;
           -2 -1 4 -1;
           -1 -2 -1 4];
% 用来拼Ma的单个矩阵
Mae = 1/36*[4 2 1 2;
            2 4 2 1;
            1 2 4 2;
            2 1 2 4];

iIndexaa = reshape(kron(edofMata,ones(4,1))',16*nelx*nely,1);
jIndexaa = reshape(kron(edofMata,ones(1,4))',16*nelx*nely,1);

disp((1/rhoa + (1/rhoav-1/rhoa)*xPhys./(1+q*(1-xPhys))));

sKa = reshape(Kae(:)*(1/rhoa + (1/rhoav-1/rhoa)*xPhys(:)./(1+q*(1-xPhys(:))))',16*nelx*nely,1);
sMa = reshape(Mae(:)*(1/kappa + (1/kappav-1/kappa)*xPhys(:))',16*nelx*nely,1);


Ka = sparse(iIndexaa,jIndexaa,sKa); Ka = (Ka+Ka')/2;
Ma = sparse(iIndexaa,jIndexaa,sMa); Ma = (Ma+Ma')/2;

disp(['Ka的行数: ', num2str(size(Ka, 1))]);  % 直接显示行数
disp(['Ka的列数: ', num2str(size(Ka, 2))]);  % 直接显示列数
disp(['Ma的行数: ', num2str(size(Ma, 1))]);  % 直接显示行数
disp(['Ma的列数: ', num2str(size(Ma, 2))]);  % 直接显示列数

% 用来拼Ke的单个矩阵
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
Kse = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
% 用来拼Me的单个矩阵
L = [4 0 2 0;
     0 4 0 2;
     2 0 4 0;
     0 2 0 4];
R = [1 0 2 0;
     0 1 0 2;
     2 0 1 0;
     0 2 0 1];
Mse = [L R;
     R' L];
Mse = Mse / 36;

iIndexss = reshape(kron(edofMats,ones(8,1))',64*nelx*nely,1);
jIndexss = reshape(kron(edofMats,ones(1,8))',64*nelx*nely,1);
sKs = reshape(Kse(:)*(Ev+(E-Ev)*xPhys(:)./(1+q*(1-xPhys(:))))',64*nelx*nely,1);
sMs = reshape(Mse(:)*(rhosv+(rhos-rhosv)*xPhys(:))',64*nelx*nely,1);

Ks = sparse(iIndexss,jIndexss,sKs); Ks = (Ks+Ks')/2;
Ms = sparse(iIndexss,jIndexss,sMs); Ms = (Ms+Ms')/2;

disp(['Ks的行数: ', num2str(size(Ks, 1))]);  % 直接显示行数
disp(['Ks的列数: ', num2str(size(Ks, 2))]);  % 直接显示列数
disp(['Ms的行数: ', num2str(size(Ms, 1))]);  % 直接显示行数
disp(['Ms的列数: ', num2str(size(Ms, 2))]);  % 直接显示列数


% 8*4
Se = [-2 0 0 -1;
       -2 -1 0 0;
        0 2 1 0;
       -1 -2 0 0;
        0 1 2 0;
        0 0 2 1;
       -1 0 0 -2;
        0 0 1 2]/6;

iIndexsa = reshape(kron(edofMats,ones(4,1))',32*nelx*nely,1);
jIndexsa = reshape(kron(edofMata,ones(1,8))',32*nelx*nely,1);

sSp = reshape(Se(:)*(1-xPhys(:)'),32*nelx*nely,1);
sSu = reshape(Se(:)*(xPhys(:)'),32*nelx*nely,1);

Sp = sparse(iIndexsa,jIndexsa,sSp);
Su = sparse(iIndexsa,jIndexsa,sSu); Su = Su';

disp(['Sp的行数: ', num2str(size(Sp, 1))]);  % 直接显示行数
disp(['Sp的列数: ', num2str(size(Sp, 2))]);  % 直接显示列数
disp(['Su的行数: ', num2str(size(Su, 1))]);  % 直接显示行数
disp(['Su的列数: ', num2str(size(Su, 2))]);  % 直接显示列数

K = [Ks Sp;
    zeros(size(Ka,1),size(Ks,2)) Ka];

M = [Ms zeros(size(Ms,1),size(Ma,2))
    Su Ma];

