%% Set Parameters
nelx=40;
nely=40;
volfrac=0.3;
penal=3.0;
rmin=1.1;

ft=2;
type=0;

%% MATERIAL PROPERTIES
E0 = 100;
Emin = 1e-9*E0;
nu = 0.3;
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
Ue = zeros(2*(nely+1)*(nelx+1),1);
Ud = zeros(2*(nely+1)*(nelx+1),1);

din=2*ceil((nely+1)/2)-1;
dout = 2*(nelx)*(nely+1)+2*ceil((nely+1)/2)-1;

L= zeros(2*(nely+1)*(nelx+1),1);    L(dout)=1;
Lambda = zeros(2*(nely+1)*(nelx+1),1);
Lambdae = zeros(2*(nely+1)*(nelx+1),1);
Lambdad = zeros(2*(nely+1)*(nelx+1),1);

F = sparse(din,1,1,2*(nely+1)*(nelx+1),1);

fixeddofs = union([1:2], [2*(nely+1)-1,2*(nely+1)]);
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
x = repmat(volfrac,nely,nelx);
xPhys = x;

beta=1;
xe = exp(-beta*(1-xPhys))-(1-xPhys)*exp(-beta);
xd = 1-exp(-beta*xPhys)+xPhys*exp(-beta);

loop = 0;
change = 1;
loopbeta = 0;
c_history = [];
ce_history = [];
cd_history = [];


if type~=0
    %% INITIALIZE MMA OPTIMIZER
    m     = 1;                % The number of general constraints. Volumn constrain
    n     = nelx*nely;        % The number of design variables x_j.
    xmin  = zeros(n,1);       % Column vector with the lower bounds for the variables x_j.
    xmax  = ones(n,1);        % Column vector with the upper bounds for the variables x_j.
    xold1 = x(:);          % xval, one iteration ago (provided that iter>1).
    xold2 = x(:);          % xval, two iterations ago (provided that iter>2).
    low   = [  ];             % Column vector with the lower asymptotes from the previous iteration (provided that iter>1).
    upp   = [  ];             % Column vector with the upper asymptotes from the previous iteration (provided that iter>1).
    a0    = 1;                % The constants a_0 in the term a_0*z.
    a     = zeros(m,1);       % Column vector with the constants a_i in the terms a_i*z.
    c_MMA = 1000*ones(m,1);   % Column vector with the constants c_i in the terms c_i*y_i.
    d     = zeros(m,1);       % Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
    move  = 0.1;              % Max allowed change in design variables
    %% START ITERATION
    while change > 0.003

        loopbeta = loopbeta+1;
        loop = loop + 1;

        sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
        K = sparse(iK,jK,sK); K = (K+K')/2;
        sKe = reshape(KE(:)*(Emin+xe(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
        Ke = sparse(iK,jK,sKe); Ke = (Ke+Ke')/2;
        sKd = reshape(KE(:)*(Emin+xd(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
        Kd = sparse(iK,jK,sKd); Kd = (Kd+Kd')/2;

        K(din,din) = K(din,din) + In_spring;
        K(dout,dout) = K(dout,dout) + Out_spring;
        Ke(din,din) = Ke(din,din) + In_spring;
        Ke(dout,dout) = Ke(dout,dout) + Out_spring;
        Kd(din,din) = Kd(din,din) + In_spring;
        Kd(dout,dout) = Kd(dout,dout) + Out_spring;

        U(freedofs) = K(freedofs,freedofs)\F(freedofs);
        Ue(freedofs) = Ke(freedofs,freedofs)\F(freedofs);
        Ud(freedofs) = Kd(freedofs,freedofs)\F(freedofs);

        Lambda(freedofs) = -K(freedofs,freedofs)\L(freedofs);
        Lambdae(freedofs) = -Ke(freedofs,freedofs)\L(freedofs);
        Lambdad(freedofs) = -Kd(freedofs,freedofs)\L(freedofs);

        %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
        c = U(dout);
        ce = Ue(dout);
        cd = Ud(dout);

        cele = reshape(sum((Lambda(edofMat)*KE).*U(edofMat),2),nely,nelx);       % ce=Lambda'*K*U
        dUdxPhys = penal*(E0-Emin)*xPhys.^(penal-1).*cele;
        celee = reshape(sum((Lambdae(edofMat)*KE).*Ue(edofMat),2),nely,nelx);
        dUedxe = penal*(E0-Emin)*xe.^(penal-1).*celee;
        celed = reshape(sum((Lambdad(edofMat)*KE).*Ud(edofMat),2),nely,nelx);
        dUddxd = penal*(E0-Emin)*xd.^(penal-1).*celed;

        dxedxPhys = beta*exp(-beta*(1-xPhys))+exp(-beta);
        dxddxPhys = beta*exp(-beta*xPhys)+exp(-beta);

        dUedxPhys = dUedxe.*dxedxPhys;
        dUddxPhys = dUddxd.*dxddxPhys;

        dv = ones(nely,nelx);

        c_history = [c_history c];
        ce_history = [ce_history ce];
        cd_history = [cd_history cd];

        %% FILTERING/MODIFICATION OF SENSITIVITIES

        dUdxPhys(:) = H*(dUdxPhys(:)./Hs);
        dUedxPhys(:) = H*(dUedxPhys(:)./Hs);
        dUddxPhys(:) = H*(dUddxPhys(:)./Hs);

        %% METHOD OF MOVING ASYMPTOTES
        xmin = max( x(:) - move, 0 );
        xmax = min( x(:) + move, 1 );

        if type == 1
            f0val = c;
            df0dx = dUdxPhys(:);
        elseif type==2
            f0val = ce;
            df0dx = dUedxPhys(:);
        elseif type==3
            f0val = cd;
            df0dx = dUddxPhys(:);
        end

        fval  = sum(xPhys(:))/(volfrac*n) - 1;
        dfdx  = dv(:)'/ (volfrac*n);

        [ xmma, ~, ~, ~, ~, ~, ~, ~, ~, low, upp ] = ...
            mmasub( m, n, loop, x(:), xmin, xmax, xold1, xold2, f0val,df0dx,fval,dfdx,low,upp,a0,a,c_MMA,d);
        % Update MMA Variables
        xold2    = xold1(:);
        xold1    = x(:);
        x     = reshape(xmma,nely,nelx);

        xPhys(:) = (H*x(:))./Hs;
        xe = exp(-beta*(1-xPhys))-(1-xPhys)*exp(-beta);
        xd = 1-exp(-beta*xPhys)+xPhys*exp(-beta);

        change = max(abs(x(:)-xold1(:)));

        %% PRINT RESULTS
        if type==1
            fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, mean(xPhys(:)),change);
            colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
        elseif type==2
            fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, mean(xe(:)),change);
            colormap(gray); imagesc(1-xe); caxis([0 1]); axis equal; axis off; drawnow;
        elseif type==3
            fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, mean(xd(:)),change);
            colormap(gray); imagesc(1-xd); caxis([0 1]); axis equal; axis off; drawnow;
        end

        if  beta < 4 && loopbeta >= 50 && change <= 0.02
            beta = 2*beta;
            loopbeta = 0;
            change = 1;
            fprintf('Parameter beta increased to %g.\n',beta);
        end
    end    % while change > 0.001

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

end

if type==0
    m     = 2+1;              % The number of general constraints. Volumn constrain
    n     = nelx*nely;        % The number of design variables x_j.
    xmin  = zeros(n,1);       % Column vector with the lower bounds for the variables x_j.
    xmax  = ones(n,1);        % Column vector with the upper bounds for the variables x_j.
    xold1 = x(:);          % xval, one iteration ago (provided that iter>1).
    xold2 = x(:);          % xval, two iterations ago (provided that iter>2).
    low   = [  ];             % Column vector with the lower asymptotes from the previous iteration (provided that iter>1).
    upp   = [  ];             % Column vector with the upper asymptotes from the previous iteration (provided that iter>1).
    a0    = 1;                % The constants a_0 in the term a_0*z.
    a     = 1*ones(m,1);      % Column vector with the constants a_i in the terms a_i*z.
    a(end)=0;
    c_MMA = 1000*ones(m,1);   % Column vector with the constants c_i in the terms c_i*y_i.
    d     = zeros(m,1);       % Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
    move  = 0.1;              % Max allowed change in design variables
    fval = zeros(m,1);
    dfdx = zeros(m,n);
    f0val = 0;
    df0dx = zeros(n, 1);

    figure;

    % 创建VideoWriter对象，用于保存视频文件
    % video = VideoWriter('structure.mp4', 'MPEG-4');  % 指定使用 MPEG-4 编码器
    % video.FrameRate = 10; % 设置帧率，例如10帧每秒
    % open(video); % 打开视频文件

    %% START ITERATION
    while change > 0.001 || beta<16

        loopbeta = loopbeta+1;
        loop = loop + 1;

        sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
        K = sparse(iK,jK,sK); K = (K+K')/2;
        sKe = reshape(KE(:)*(Emin+xe(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
        Ke = sparse(iK,jK,sKe); Ke = (Ke+Ke')/2;
        sKd = reshape(KE(:)*(Emin+xd(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
        Kd = sparse(iK,jK,sKd); Kd = (Kd+Kd')/2;

        K(din,din) = K(din,din) + In_spring;
        K(dout,dout) = K(dout,dout) + Out_spring;
        Ke(din,din) = Ke(din,din) + In_spring;
        Ke(dout,dout) = Ke(dout,dout) + Out_spring;
        Kd(din,din) = Kd(din,din) + In_spring;
        Kd(dout,dout) = Kd(dout,dout) + Out_spring;

        U(freedofs) = K(freedofs,freedofs)\F(freedofs);
        Ue(freedofs) = Ke(freedofs,freedofs)\F(freedofs);
        Ud(freedofs) = Kd(freedofs,freedofs)\F(freedofs);

        Lambda(freedofs) = -K(freedofs,freedofs)\L(freedofs);
        Lambdae(freedofs) = -Ke(freedofs,freedofs)\L(freedofs);
        Lambdad(freedofs) = -Kd(freedofs,freedofs)\L(freedofs);

        %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
        c = U(dout);
        ce = Ue(dout);
        cd = Ud(dout);

        cele = reshape(sum((Lambda(edofMat)*KE).*U(edofMat),2),nely,nelx);       % ce=Lambda'*K*U
        dUdxPhys = penal*(E0-Emin)*xPhys.^(penal-1).*cele;
        celee = reshape(sum((Lambdae(edofMat)*KE).*Ue(edofMat),2),nely,nelx);
        dUedxe = penal*(E0-Emin)*xe.^(penal-1).*celee;
        celed = reshape(sum((Lambdad(edofMat)*KE).*Ud(edofMat),2),nely,nelx);
        dUddxd = penal*(E0-Emin)*xd.^(penal-1).*celed;

        dxedxPhys = beta*exp(-beta*(1-xPhys))+exp(-beta);
        dxddxPhys = beta*exp(-beta*xPhys)+exp(-beta);

        dUedxPhys = dUedxe.*dxedxPhys;
        dUddxPhys = dUddxd.*dxddxPhys;

        dv = ones(nely,nelx);

        c_history = [c_history c];
        ce_history = [ce_history ce];
        cd_history = [cd_history cd];

        %% FILTERING/MODIFICATION OF SENSITIVITIES

        dUdxPhys(:) = H*(dUdxPhys(:)./Hs);
        dUedxPhys(:) = H*(dUedxPhys(:)./Hs);
        dUddxPhys(:) = H*(dUddxPhys(:)./Hs);

        %% METHOD OF MOVING ASYMPTOTES
        xmin = max( x(:) - move, 0 );
        xmax = min( x(:) + move, 1 );

        fval(1) = 10+c;
        fval(2) = 10+ce;
        fval(3) = 10+cd;
        fval(m) = sum(xPhys(:))/(volfrac*n) - 1;

        dfdx(1,:) = reshape(dUdxPhys,1,n);
        dfdx(2,:) = reshape(dUedxPhys,1,n);
        dfdx(3,:) = reshape(dUddxPhys,1,n);
        dfdx(m,:) = dv(:)'/ (volfrac*n);

        [ xmma, ~, ~, ~, ~, ~, ~, ~, ~, low, upp ] = ...
            mmasub( m, n, loop, x(:), xmin, xmax, xold1, xold2, f0val,df0dx,fval,dfdx,low,upp,a0,a,c_MMA,d);
        % Update MMA Variables
        xold2    = xold1(:);
        xold1    = x(:);
        x     = reshape(xmma,nely,nelx);

        xPhys(:) = (H*x(:))./Hs;
        xe = exp(-beta*(1-xPhys))-(1-xPhys)*exp(-beta);
        xd = 1-exp(-beta*xPhys)+xPhys*exp(-beta);

        change = max(abs(x(:)-xold1(:)));

        %% PRINT RESULTS
        fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, mean(xPhys(:)),change);


        % 假设已经创建了 figure，继续添加子图
        % 第一个图：1 - xPhys
        subplot('Position', [0.05, 0.1, 0.27, 0.8]);  % 调整 [左, 下, 宽, 高]
        colormap(gray);
        imagesc(1 - xPhys);
        caxis([0 1]);
        axis equal;
        axis off;
        title(['original, c = ', num2str(c, '%.4f')], 'FontSize', 12);  % 确保标题在图片上方正中

        % 第二个图：1 - xe
        subplot('Position', [0.37, 0.1, 0.27, 0.8]);  % 调整 [左, 下, 宽, 高]
        colormap(gray);
        imagesc(1 - xe);
        caxis([0 1]);
        axis equal;
        axis off;
        title(['erode, ce = ', num2str(ce, '%.4f')], 'FontSize', 12);  % 确保标题在图片上方正中

        % 第三个图：1 - xd
        subplot('Position', [0.69, 0.1, 0.27, 0.8]);  % 调整 [左, 下, 宽, 高]
        colormap(gray);
        imagesc(1 - xd);
        caxis([0 1]);
        axis equal;
        axis off;
        title(['dilate, cd = ', num2str(cd, '%.4f')], 'FontSize', 12);  % 确保标题在图片上方正中

        drawnow;  % 更新绘图窗口

        % 将当前帧写入视频
        % frame = getframe(gcf);  % 获取当前figure的帧
        % writeVideo(video, frame);  % 将帧写入视频文件

        if  beta < 16 && loopbeta >= 50 && change <= 0.02
            beta = 2*beta;
            loopbeta = 0;
            change = 1;
            fprintf('Parameter beta increased to %g.\n',beta);
        end
    end    % while change > 0.001

    saveas(gcf, 'structure.png');

    % 关闭视频文件
    % close(video);

    figure('Position', [100, 100, 800, 600]);  % 调整图像的大小，设置宽度为800，高度为600

    % 第四张图：绘制 c_history、ce_history 和 cd_history 在同一张图上
    subplot('Position', [0.1, 0.1, 0.8, 0.8]);  % 调整大小和位置 [左, 下, 宽, 高]
    plot(1:loop, c_history, '-o', 'LineWidth', 2);  % 绘制 c_history
    hold on;
    plot(1:loop, ce_history, '-x', 'LineWidth', 2);  % 绘制 ce_history
    plot(1:loop, cd_history, '-s', 'LineWidth', 2);  % 绘制 cd_history
    hold off;
    xlabel('Iteration');
    ylabel('Objective Function Value');
    legend('c', 'ce', 'cd');  % 添加图例
    title('Objective Function History (Multiple)');

    save_filename = sprintf('displacement.png');
    saveas(gcf, save_filename);
end % if type==0



