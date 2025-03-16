%% 参数设置
c = 340.27;
kappa = 141834.999; % 气体的体积模量
kappav = kappa*10^7; % 搭配气体的虚拟的体积模量
rhoa = 1.225; % 气体的密度，待确定
rhoav = rhoa*10^6; % 搭配气体的虚拟的密度
k = 1.5239; % 似乎是波数，等于omega/c_a
omega = 518.55;
pin = 10;
Lx = 2;
Ly = 1;



nelx = 400;  % x方向网格数
nely = 200;  % y方向网格数

edge = 100/nely;







xPhys = zeros(nely, nelx);

q = 1; % RAMP插值的参数



% 声压场的编号，整个区域的声压场都需要计算，一共(nely+1)*(nelx+1)个节点，也是这么多自由度
P = zeros((nely+1)*(nelx+1),1); % 用于存储声压场的自由度
nodenrsa = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx); % 给(nely+1)*(nelx+1)个节点编号
edofVeca = reshape(nodenrsa(1:end-1,1:end-1)+nely+1,nelx*nely,1); % 记录每个单元的右上节点的编号
edofMata = repmat(edofVeca,1,4)+repmat([0 -nely-1 -nely 1],nelx*nely,1); % 记录每个单元的四个声压自由度的编号
% 以上这部分在纯声学计算特征频率的代码中用过，所以应该是没错的


%% PREPARE FINITE ELEMENT ANALYSIS
% 纯声学部分
Kae = 1/6*[4 -1 -2 -1; -1 4 -1 -2; -2 -1 4 -1; -1 -2 -1 4];
Mae = 1/36*[4 2 1 2; 2 4 2 1; 1 2 4 2; 2 1 2 4];

iIndexaa = reshape(kron(edofMata,ones(4,1))',16*nelx*nely,1);
jIndexaa = reshape(kron(edofMata,ones(1,4))',16*nelx*nely,1);

Mle = 1/6*[2 0 0 1;0 0 0 0;0 0 0 0;1 0 0 2];
Mre = 1/6*[0 0 0 0;0 2 1 0;0 1 2 0;0 0 0 0];
iIndexl = reshape(kron(edofMata(1:nely, :),ones(4,1))',16*nely,1);
jIndexl = reshape(kron(edofMata(1:nely, :),ones(1,4))',16*nely,1);
iIndexr = reshape(kron(edofMata((nelx-1)*nely+1:nelx*nely, :),ones(4,1))',16*nely,1);
jIndexr = reshape(kron(edofMata((nelx-1)*nely+1:nelx*nely, :),ones(1,4))',16*nely,1);
% 这部分比较简单，感觉不会错

% 纯声学部分
sKa = reshape(Kae(:)*(1/rhoa + (1/rhoav-1/rhoa)*xPhys(:)./(1+q*(1-xPhys(:))))',16*nelx*nely,1);
sMa = reshape(Mae(:)*(1/kappa + (1/kappav-1/kappa)*xPhys(:))',16*nelx*nely,1);
% 用的是论文中的插值公式，除非我对matlab的语法理解有问题，不然不会错
Ka = sparse(iIndexaa,jIndexaa,sKa); Ka = (Ka+Ka')/2;
Ma = sparse(iIndexaa,jIndexaa,sMa); Ma = (Ma+Ma')/2;

sMl = reshape(Mle(:) * (1/rhoa + (1/rhoav - 1/rhoa) * xPhys(:,1) ./ (1 + q*(1 - xPhys(:,1))))', 16*nely, 1);
sMr = reshape(Mre(:) * (1/rhoa + (1/rhoav - 1/rhoa) * xPhys(:,nelx) ./ (1 + q*(1 - xPhys(:,nelx))))', 16*nely, 1);
% sMl = reshape(Mle(:)*(1/rhoa + (1/rhoav-1/rhoa)*repmat(xPhys(:,1), size(xPhys,2),1)./(1+q*(1-repmat(xPhys(:,1), size(xPhys,2),1))))',16*nely,1);
% sMr = reshape(Mre(:)*(1/rhoa + (1/rhoav-1/rhoa)*repmat(xPhys(:,nelx), size(xPhys,2),1)./(1+q*(1-repmat(xPhys(:,nelx), size(xPhys,2),1))))',16*nely,1);
Ml = sparse(iIndexl,jIndexl,sMl,(nely+1)*(nelx+1),(nely+1)*(nelx+1)); Ml = (Ml+Ml')/2;
Mr = sparse(iIndexr,jIndexr,sMr); Mr = (Mr+Mr')/2;
Mi = Ml+Mr;

Fa = zeros((nely+1)*(nelx+1), 1);
Fa(1:nely+1) = [(1/rhoa + (1/rhoav-1/rhoa)*xPhys(:,1)./(1+q*(1-xPhys(:,1)))); 0]/2 + [0; (1/rhoa + (1/rhoav-1/rhoa)*xPhys(:,1)./(1+q*(1-xPhys(:,1))))]/2;

K = edge^2*Ka-edge^2*omega^2*Ma+edge*1i*k*Ml;
F = edge*2*1i*k*pin*Fa;
P = reshape(K\F,nely+1,nelx+1);

%% 提取中间行并绘图
    mid_row = round((nely+1)/2);       % 计算中间行索引
    x_coords = linspace(0, nelx, nelx+1); % 生成x坐标
    plot(x_coords, real(P(mid_row,:)), 'LineWidth',1.5,...
        'DisplayName',sprintf('nely=%d, nelx=%d',nely,nelx));