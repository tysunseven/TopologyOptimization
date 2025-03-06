%% 参数设置
alpha = 7;
E = 1000; % 固体的杨氏模量
Ev = E/10^(alpha-1); % 搭配固体的虚拟的杨氏模量
nu = 0.3; % 固体的泊松比
rhos = 1; % 固体的密度
rhosv = rhos/10^alpha; % 搭配固体的虚拟的密度

kappa = 1; % 气体的体积模量
kappav = kappa*10^alpha; % 搭配气体的虚拟的体积模量
rhoa = 1; % 气体的密度，待确定
rhoav = rhoa*10^(alpha-1); % 搭配气体的虚拟的密度

q = 1; % RAMP插值的参数
omega = 2;

%% 设计变量的初始化
nelxs = 25;
nelx = nelxs*9;  % x方向网格数,4.5米,每个单元格0.02米
nely = nelxs*4;  % y方向网格数,2米,每个单元格0.02米
volfrac = 0.55; % 体积分数，也决定了初始解

x = zeros(nely, nelx);   % 初始化全0矩阵
x(:,nelxs*4+1:nelxs*5) = volfrac;    % 将中间25列的值设为0.55
xPhys = x;  % 过滤和Heaviside投影在之后补充
xPhyss = xPhys(:,nelxs*4+1:nelxs*5); %把中间部分取出来
% imagesc(1-x); colormap(gray); clim([0 1]); axis equal; axis off;  % 测试代码

%% 物理场的编号
% 声压场的编号，整个区域的声压场都需要计算，一共(nely+1)*(nelx+1)个节点，也是这么多自由度
P = zeros((nely+1)*(nelx+1),1); % 用于存储声压场的自由度
nodenrsa = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx); % 给(nely+1)*(nelx+1)个节点编号
edofVeca = reshape(nodenrsa(1:end-1,1:end-1)+nely+1,nelx*nely,1); % 记录每个单元的右上节点的编号
edofMata = repmat(edofVeca,1,4)+repmat([0 -nely-1 -nely 1],nelx*nely,1); % 记录每个单元的四个声压自由度的编号
% 以上这部分在纯声学计算特征频率的代码中用过，所以应该是没错的

% 位移场的编号，中间区域的位移场需要计算，一共(nely+1)*(nelxs+1)个节点，二倍的自由度
U = zeros(2*(nely+1)*(nelxs+1),1); % 用于存储位移场的自由度
nodenrss = reshape(1:(1+nelxs)*(1+nely),1+nely,1+nelxs); % 给(nely+1)*(nelx+1)个节点编号
edofVecs = reshape(2*nodenrss(1:end-1,1:end-1)+1,nelxs*nely,1); % 记录每个单元的左下节点的编号
edofMats = repmat(edofVecs,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelxs*nely,1); % 记录每个单元的八个位移自由度的编号
% 以上这部分在纯力学的代码中用过很多次了，所以应该是没错的

%% 局部矩阵的计算
% 纯声学部分
Kae = 1/6*[4 -1 -2 -1; -1 4 -1 -2; -2 -1 4 -1; -1 -2 -1 4];
Mae = 1/36*[4 2 1 2; 2 4 2 1; 1 2 4 2; 2 1 2 4];

iIndexaa = reshape(kron(edofMata,ones(4,1))',16*nelx*nely,1);
jIndexaa = reshape(kron(edofMata,ones(1,4))',16*nelx*nely,1);
% 以上这部分在纯声学计算特征频率的代码中用过，所以应该是没错的

% 纯力学部分
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
Kse = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);

L = [4 0 2 0; 0 4 0 2; 2 0 4 0; 0 2 0 4];
R = [1 0 2 0; 0 1 0 2; 2 0 1 0; 0 2 0 1];
Mse = [L R; R' L]/36;

iIndexss = reshape(kron(edofMats,ones(8,1))',64*nelxs*nely,1);
jIndexss = reshape(kron(edofMats,ones(1,8))',64*nelxs*nely,1);
% 唯一需要检查的就是Mse是不是编对了，但他是在Mae基础上来的，大概率对


% 声振耦合的部分
Spe = [-2 0 0 -1; -2 -1 0 0; 0 2 1 0; -1 -2 0 0; 0 1 2 0; 0 0 2 1; -1 0 0 -2; 0 0 1 2]/6;
Sue = Spe';
% 这部分仔细检查验算了，应该是对的

iIndexsa = reshape(kron(edofMats,ones(4,1))',32*nelxs*nely,1);
jIndexsa = reshape(kron(edofMata,ones(1,8))',32*nelxs*nely,1);
% 感觉应该是要把edofMata改成中间需要的一些行

%% 全局矩阵的拼装
% 纯声学部分
sKa = reshape(Kae(:)*(1/rhoa + (1/rhoav-1/rhoa)*xPhys(:)./(1+q*(1-xPhys(:))))',16*nelx*nely,1);
sMa = reshape(Mae(:)*(1/kappa + (1/kappav-1/kappa)*xPhys(:))',16*nelx*nely,1);
% 用的是论文中的插值公式，除非我对matlab的语法理解有问题，不然不会错
Ka = sparse(iIndexaa,jIndexaa,sKa); Ka = (Ka+Ka')/2;
Ma = sparse(iIndexaa,jIndexaa,sMa); Ma = (Ma+Ma')/2;

% 纯力学部分
sKs = reshape(Kse(:)*(Ev+(E-Ev)*xPhyss(:)./(1+q*(1-xPhyss(:))))',64*nelxs*nely,1);
sMs = reshape(Mse(:)*(rhosv+(rhos-rhosv)*xPhyss(:))',64*nelxs*nely,1);
% 用的是论文中的插值公式，除非我对matlab的语法理解有问题，不然不会错
Ks = sparse(iIndexss,jIndexss,sKs); Ks = (Ks+Ks')/2;
Ms = sparse(iIndexss,jIndexss,sMs); Ms = (Ms+Ms')/2;

% 声振耦合部分
sSp = reshape(Spe(:)*(xPhyss(:)'),32*nelxs*nely,1);
sSu = -reshape(Sue(:)*(xPhyss(:)'),32*nelxs*nely,1);

Sp = sparse(iIndexsa,jIndexsa,sSp);
Su = sparse(iIndexsa,jIndexsa,sSu);

K = [Ks Sp;
    zeros(size(Ka,1),size(Ks,2)) Ka];

M = [Ms zeros(size(Ms,1),size(Ma,2))
    Su Ma];

%% 有限元的计算




