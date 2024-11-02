% Function: Multi-cell microstructure optimization with compatibility
% Author: Eric Garner (e.garner@tudelft.nl)
% Version: 2018-11-12

% Example  topX_Dual(3,100,0.2,0.6,0,1,1,0,0,1,1000)

function [xPhys_cat, c, Q] = topX_Dual(CellCount,nel,Vmin,Vmax,w1,OBJECTIVE,GLOBAL,LOCAL,ISO,FilterDom,nloop)


%% INPUTS
% OBJECTIVE 1: Bulk Modulus; 2: Shear Modulus; 3: Negative Poisson
% CONSTRAINT TOGGLES (GLOBAL, LOCAL, ISO): enabled(1), disabled(0)
% GLOBAL SCALE FILTERING TOGGLE (FilterDom): enabled(1), disabled(0)

%% ROBUST FORMULATION PARAMETERS
eta = 0.3 ;
beta = 1;
BetaUpdFreq = 200;
BetaUpdFactor = 2;

%% ISO CONSTRAINT PARAMETERS
epsilon_iso = 0.001;
EpsilonUpdFreq = 500;
EpsilonUpdFactor = 0.5;

%% LOCAL VOLUME CONSTRAINRT PARAMETERS
p = 16;                         % pNorm
r_hat = 10;                     % pNorm radius: sets maximum length scale (solid material)
v = 0.8;

%% MATERIAL PARAMETERS
E0 = 1;
nu = 0.3;
Emin = 1e-9;

%% MISC PARAMETERS
penal = 5;
rmin = 2.2;
sym = 4;                                % 0: None 1: 4-fold diagonal 2: 4-fold normal 3: Chiral 4: Cubic
Influence = CellCount;                  % Order of influence; 1: Individual, 2: Adjacent cell compatibility
IC = 1;

if w1 == 1
    Influence = 1;
end

nelx = nel;
nely = nel;

CaseCount = Influence*CellCount - Influence*(Influence-1)/2;
volfrac = linspace(Vmin,Vmax,CellCount);
vol_max = linspace(v,v,CellCount);
vol_max_pNorm = (nelx*nely*vol_max.^p).^(1/p);
w = [w1,1-w1];
w = 100*w;

% used for importing existing design
x_imp = 0;
imp = 0;

%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);

nodenrs = cell(2,1);
edofMat = cell(2,1);
iK = cell(2,1);
jK = cell(2,1);
[nodenrs{1}, edofMat{1}, iK{1}, jK{1}] = PrepareFEA(nelx,nely);
[nodenrs{2}, edofMat{2}, iK{2}, jK{2}] = PrepareFEA(2*nelx,nely);


%% PREPARE PDE FILTER
edofVecF = reshape(nodenrs{1}(1:end-1,1:end-1),nelx*nely,1);
edofMatF = repmat(edofVecF,1,4)+repmat([0 nely+[1:2] 1],nelx*nely,1);
iKF = reshape(kron(edofMatF,ones(4,1))',16*nelx*nely,1);
jKF = reshape(kron(edofMatF,ones(1,4))',16*nelx*nely,1);
iTF = reshape(edofMatF,4*nelx*nely,1);
jTF = reshape(repmat([1:nelx*nely],4,1)',4*nelx*nely,1);
sTF = repmat(1/4,4*nelx*nely,1);
TF = sparse(iTF,jTF,sTF);

Rmin = r_hat/2/sqrt(3);
KEF = Rmin^2*[4 -1 -2 -1; -1  4 -1 -2; -2 -1  4 -1; -1 -2 -1  4]/6 + ...
    [4  2  1  2;  2  4  2  1;  1  2  4  2;  2  1  2  4]/36;
sKF = reshape(KEF(:)*ones(1,nelx*nely),16*nelx*nely,1);
KF = sparse(iKF,jKF,sKF);
LF = chol(KF,'lower');

%% PERIODIC BOUNDARY CONDITIONS
U_e = cell(2,1);
U_i = cell(2,1);
U_d = cell(2,1);
ufixed = cell(2,1);
d1 = cell(2,1);
d2 = cell(2,1);
d3 = cell(2,1);
d4 = cell(2,1);
wfixed = cell(2,1);
[U_e{1}, ufixed{1}, d1{1}, d2{1}, d3{1}, d4{1}, wfixed{1}] = ApplyPeriodicBC(nelx,nely,nodenrs{1,1});
[U_i{1}, ufixed{1}, d1{1}, d2{1}, d3{1}, d4{1}, wfixed{1}] = ApplyPeriodicBC(nelx,nely,nodenrs{1,1});
[U_d{1}, ufixed{1}, d1{1}, d2{1}, d3{1}, d4{1}, wfixed{1}] = ApplyPeriodicBC(nelx,nely,nodenrs{1,1});
[U_e{2}, ufixed{2}, d1{2}, d2{2}, d3{2}, d4{2}, wfixed{2}] = ApplyPeriodicBC(2*nelx,nely,nodenrs{2,1});
[U_i{2}, ufixed{2}, d1{2}, d2{2}, d3{2}, d4{2}, wfixed{2}] = ApplyPeriodicBC(2*nelx,nely,nodenrs{2,1});
[U_d{2}, ufixed{2}, d1{2}, d2{2}, d3{2}, d4{2}, wfixed{2}] = ApplyPeriodicBC(2*nelx,nely,nodenrs{2,1});

%% INITIALIZE ITERATION
Q_e = cell(CaseCount,1);
Q_i = cell(CaseCount,1);
Q_d = cell(CaseCount,1);
dQ_e = cell(CaseCount,1);
dQ_i = cell(CaseCount,1);
dQ_d = cell(CaseCount,1);

x = cell(CaseCount,1);
xTilde = cell(CaseCount,1);
xPhys_e = cell(CaseCount,1);
xPhys_i = cell(CaseCount,1);
xPhys_d = cell(CaseCount,1);
xold1 = cell(CaseCount,1);
xold2 = cell(CaseCount,1);

for i = 1:CaseCount
    Q_e{i,1} = zeros(3,3);
    Q_i{i,1} = zeros(3,3);
    Q_d{i,1} = zeros(3,3);
    dQ_e{i,1} = cell(3,3);
    dQ_i{i,1} = cell(3,3);
    dQ_d{i,1} = cell(3,3);
end

for i = 1:CellCount
    x_imp_temp = zeros(nely,nelx);
    if imp == 1
        x_imp_temp = x_imp(:,(1+(i-1)*nelx:i*nelx));
    end
    [x{i}, xTilde{i}, xPhys_e{i}, xPhys_i{i}, xPhys_d{i}, xold1{i}, xold2{i}] = InitializeIteration(nelx,nely,beta,eta,volfrac(i),IC,x_imp_temp,imp);
end

i = CellCount + 1;
l = 1;
while i <= CaseCount
    j = 1;
    k = j + l;
    while k <= CellCount
        xTilde{i} = horzcat(xTilde{j},xTilde{k});
        xPhys_e{i} = horzcat(xPhys_e{j},xPhys_e{k});
        xPhys_i{i} = horzcat(xPhys_i{j},xPhys_i{k});
        xPhys_d{i} = horzcat(xPhys_d{j},xPhys_d{k});
        xold1{i} = horzcat(xold1{j},xold1{k});
        xold2{i} = horzcat(xold2{j},xold2{k});
        
        i  = i + 1;
        j = j + 1;
        k = j + l;
    end
    l = l + 1;
end

loop = 0;
loopbeta = 0;
change = 1;
low = 0;
upp = 0;


%% START ITERATION
% store results
c_hist = zeros(nloop,CaseCount);        % objective
vol_hist = zeros(nloop,CaseCount);      % volume
change_hist = zeros(nloop,1);   % maximum design change
sharp_hist = zeros(nloop,1);    % sharpness
cons_hist = zeros(nloop,2);     % constraints

%% START ITERATION
while (change > 0.001) && loop<nloop
    loop = loop + 1;
    loopbeta = loopbeta + 1;
    
    %% PREPARE FILTER
    iH = ones(CellCount*nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
    jH = ones(size(iH));
    sH = zeros(size(iH));
    k = 0;
    for i1 = 1:CellCount*nelx
        for j1 = 1:nely
            e1 = (i1-1)*nely+j1;
            for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),CellCount*nelx)
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
    H_ind = sparse(iH,jH,sH);
    Hs_ind = sum(H_ind,2);
    
    %% FE ANALYSIS
    for i = 1:CaseCount
        if i <= CellCount
            U_e{i} = FEAnalysis(nelx,nely,KE,Emin,xPhys_e{i},penal,E0,iK{1},jK{1},d1{1},d2{1},d3{1},d4{1},ufixed{1},wfixed{1});
            U_i{i} = FEAnalysis(nelx,nely,KE,Emin,xPhys_i{i},penal,E0,iK{1},jK{1},d1{1},d2{1},d3{1},d4{1},ufixed{1},wfixed{1});
            U_d{i} = FEAnalysis(nelx,nely,KE,Emin,xPhys_d{i},penal,E0,iK{1},jK{1},d1{1},d2{1},d3{1},d4{1},ufixed{1},wfixed{1});
            [Q_e{i}, dQ_e{i}] = Homogenization(nelx,nely,U_e{i},edofMat{1},Emin,E0,xPhys_e{i},penal,KE);
            [Q_i{i}, dQ_i{i}] = Homogenization(nelx,nely,U_i{i},edofMat{1},Emin,E0,xPhys_i{i},penal,KE);
            [Q_d{i}, dQ_d{i}] = Homogenization(nelx,nely,U_d{i},edofMat{1},Emin,E0,xPhys_d{i},penal,KE);
        else
            U_e{i} = FEAnalysis(2*nelx,nely,KE,Emin,xPhys_e{i},penal,E0,iK{2},jK{2},d1{2},d2{2},d3{2},d4{2},ufixed{2},wfixed{2});
            U_i{i} = FEAnalysis(2*nelx,nely,KE,Emin,xPhys_i{i},penal,E0,iK{2},jK{2},d1{2},d2{2},d3{2},d4{2},ufixed{2},wfixed{2});
            U_d{i} = FEAnalysis(2*nelx,nely,KE,Emin,xPhys_d{i},penal,E0,iK{2},jK{2},d1{2},d2{2},d3{2},d4{2},ufixed{2},wfixed{2});
            [Q_e{i}, dQ_e{i}] = Homogenization(2*nelx,nely,U_e{i},edofMat{2},Emin,E0,xPhys_e{i},penal,KE);
            [Q_i{i}, dQ_i{i}] = Homogenization(2*nelx,nely,U_i{i},edofMat{2},Emin,E0,xPhys_i{i},penal,KE);
            [Q_d{i}, dQ_d{i}] = Homogenization(2*nelx,nely,U_d{i},edofMat{2},Emin,E0,xPhys_d{i},penal,KE);
        end
    end
    
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    c_e = 0;
    c_i = 0;
    c_d = 0;
    dc_e = zeros(nely,CellCount*nelx);
    dc_i = zeros(nely,CellCount*nelx);
    dc_d = zeros(nely,CellCount*nelx);
    
    for i = 1:CellCount
        
        [c_e_temp, dc_e_temp] = Objective(Q_e{i},dQ_e{i},OBJECTIVE,loop);
        [c_i_temp,dc_i_temp] = Objective(Q_i{i},dQ_i{i},OBJECTIVE,loop);
        [c_d_temp,dc_d_temp] = Objective(Q_d{i},dQ_d{i},OBJECTIVE,loop);
        
        c_hist(loop,i) = c_i_temp;
        
        if i == 1 || i == CellCount
            c_e =  c_e + 0.5 * w(1) * c_e_temp;
            dc_e(1:nely,(((i-1)*nelx)+1):(i*nelx)) = dc_e(1:nely,(((i-1)*nelx)+1):(i*nelx)) + 0.5*w(1) * dc_e_temp;
            c_i =  c_i + 0.5 * w(1) * c_i_temp;
            dc_i(1:nely,(((i-1)*nelx)+1):(i*nelx)) = dc_i(1:nely,(((i-1)*nelx)+1):(i*nelx)) + 0.5*w(1) * dc_i_temp;
            c_d =  c_d + 0.5 * w(1) * c_d_temp;
            dc_d(1:nely,(((i-1)*nelx)+1):(i*nelx)) = dc_d(1:nely,(((i-1)*nelx)+1):(i*nelx)) + 0.5*w(1) * dc_d_temp;
        end
        c_e =  c_e + w(1) * c_e_temp;
        dc_e(1:nely,(((i-1)*nelx)+1):(i*nelx)) = dc_e(1:nely,(((i-1)*nelx)+1):(i*nelx)) + w(1) * dc_e_temp;
        c_i =  c_i + w(1) * c_i_temp;
        dc_i(1:nely,(((i-1)*nelx)+1):(i*nelx)) = dc_i(1:nely,(((i-1)*nelx)+1):(i*nelx)) + w(1) * dc_i_temp;
        c_d =  c_d + w(1) * c_d_temp;
        dc_d(1:nely,(((i-1)*nelx)+1):(i*nelx)) = dc_d(1:nely,(((i-1)*nelx)+1):(i*nelx)) + w(1) * dc_d_temp;
    end
    
    i = CellCount + 1;
    l = 1;
    while i <= CaseCount
        j = 1;
        k = j + l;
        while k <= CellCount
            
            [c_e_temp, dc_e_temp] = Objective(Q_e{i},dQ_e{i},OBJECTIVE,loop);
            [c_i_temp,dc_i_temp] = Objective(Q_i{i},dQ_i{i},OBJECTIVE,loop);
            [c_d_temp,dc_d_temp] = Objective(Q_d{i},dQ_d{i},OBJECTIVE,loop);
            
            c_hist(loop,i) = c_i_temp;
            
            c_e = c_e + w(2) * c_e_temp;
            tmp_e = + w(2) * dc_e_temp;
            tmp_e_L = tmp_e(1:nely,1:nelx);
            tmp_e_R = tmp_e(1:nely,(nelx+1):end);
            dc_e(1:nely,(((j-1)*nelx)+1):(j*nelx)) = dc_e(1:nely,(((j-1)*nelx)+1):(j*nelx)) + tmp_e_L;
            dc_e(1:nely,(((k-1)*nelx)+1):(k*nelx)) = dc_e(1:nely,(((k-1)*nelx)+1):(k*nelx)) + tmp_e_R;
            
            c_i = c_i + w(2) * c_i_temp;
            tmp_i = + w(2) * dc_i_temp;
            tmp_i_L = tmp_i(1:nely,1:nelx);
            tmp_i_R = tmp_i(1:nely,(nelx+1):end);
            dc_i(1:nely,(((j-1)*nelx)+1):(j*nelx)) = dc_i(1:nely,(((j-1)*nelx)+1):(j*nelx)) + tmp_i_L;
            dc_i(1:nely,(((k-1)*nelx)+1):(k*nelx)) = dc_i(1:nely,(((k-1)*nelx)+1):(k*nelx)) + tmp_i_R;
            
            c_d = c_d + w(2) * c_d_temp;
            tmp_d = + w(2) * dc_d_temp;
            tmp_d_L = tmp_d(1:nely,1:nelx);
            tmp_d_R = tmp_d(1:nely,(nelx+1):end);
            dc_d(1:nely,(((j-1)*nelx)+1):(j*nelx)) = dc_d(1:nely,(((j-1)*nelx)+1):(j*nelx)) + tmp_d_L;
            dc_d(1:nely,(((k-1)*nelx)+1):(k*nelx)) = dc_d(1:nely,(((k-1)*nelx)+1):(k*nelx)) + tmp_d_R;
            
            i  = i + 1;
            j = j + 1;
            k = j + l;
        end
        l = l + 1;
    end
    
    dv = ones(nely,CellCount*nelx);
    
    x_cat = horzcat(x{1:CellCount});
    xTilde_cat = horzcat(xTilde{1:CellCount});
    xPhys_d_cat = horzcat(xPhys_d{1:CellCount});
    
    x_pde_hat_d = cell(CellCount,1);
    dfdx_pde_d = cell(CellCount,1);
    for i = 1:CellCount
        x_pde_hat_d{i} = (TF'*(LF'\(LF\(TF*xPhys_d{i}(:)))));
        dfdx_pde_d{i} = (sum(x_pde_hat_d{i}.^p))^(1/p-1) * x_pde_hat_d{i}.^(p-1);
    end
    
    %% FILTERING AND MODIFICATION OF SENSITIVITIES
    dx_ind = cell(CellCount,1);
    for i = 1:CellCount
        dx_ind{i} = beta * (1-tanh(beta*(xTilde{i}-0.5)).*tanh(beta*(xTilde{i}-0.5))) / (tanh(beta*0.5) + tanh(beta*(1-0.5)));
    end
    
    if FilterDom == 1
        dx = beta * (1-tanh(beta*(xTilde_cat-0.5)).*tanh(beta*(xTilde_cat-0.5))) / (tanh(beta*0.5) + tanh(beta*(1-0.5)));
        dc_e(:) = H*(dc_e(:).*dx(:)./Hs);
        dc_i(:) = H*(dc_i(:).*dx(:)./Hs);
        dc_d(:) = H*(dc_d(:).*dx(:)./Hs);
        dv(:) = H*(dv(:).*dx(:)./Hs);
    else
        for i = 1:CellCount
            dc_e_temp = dc_e(1:nely,(((i-1)*nelx)+1):(i*nelx));
            dc_e_temp(:) = H_ind*(dc_e_temp(:).*dx_ind{i}(:)./Hs_ind);
            dc_e(1:nely,(((i-1)*nelx)+1):(i*nelx)) = dc_e_temp;
            
            dc_i_temp = dc_i(1:nely,(((i-1)*nelx)+1):(i*nelx));
            dc_i_temp(:) = H_ind*(dc_i_temp(:).*dx_ind{i}(:)./Hs_ind);
            dc_i(1:nely,(((i-1)*nelx)+1):(i*nelx)) = dc_i_temp;
            
            dc_d_temp = dc_d(1:nely,(((i-1)*nelx)+1):(i*nelx));
            dc_d_temp(:) = H_ind*(dc_d_temp(:).*dx_ind{i}(:)./Hs_ind);
            dc_d(1:nely,(((i-1)*nelx)+1):(i*nelx)) = dc_d_temp;
            
            dv_temp = dv(1:nely,(((i-1)*nelx)+1):(i*nelx));
            dv_temp(:) = H_ind*(dv_temp(:).*dx_ind{i}(:)./Hs_ind);
            dv(1:nely,(((i-1)*nelx)+1):(i*nelx)) = dv_temp;
        end
    end
    
    %% MMA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
    
    m = (GLOBAL + LOCAL + ISO)*CellCount;
    n = CellCount*nelx*nely;
    
    move = 0.1;
    fval = zeros((GLOBAL+LOCAL+ISO)*CellCount,1);
    dfdx = zeros((GLOBAL+LOCAL+ISO)*CellCount,nelx*nely);
    
    c = max([c_e,c_i,c_d]);
    f0val = c;
    if f0val == c_e
        df0dx = reshape(dc_e,[CellCount*nelx*nely,1]);
        Q = Q_e;
        dQ = dQ_e;
        
    elseif f0val == c_i
        df0dx = reshape(dc_i,[CellCount*nelx*nely,1]);
        Q = Q_i;
        dQ = dQ_i;
        
    elseif f0val == c_d
        df0dx = reshape(dc_d,[CellCount*nelx*nely,1]);
        Q = Q_d;
        dQ = dQ_d;
    end
    
    df0dx2 = 0*df0dx;
    
    if mod(loop,20) == 1
        volfrac_d = zeros(CellCount,1);
        for i = 1:CellCount
            volfrac_d(i) = volfrac(i)*sum(sum(xPhys_d{i}))/sum(sum(xPhys_i{i}));
        end
    end
    
    M = 0;
    if GLOBAL == 1
        for i = 1:CellCount
            fval(M*CellCount+i) = 1000*(sum(sum(xPhys_d_cat(1:nely,1+((i-1)*nelx):i*nelx))) / (nelx*nely*volfrac_d(i)) - 1);
            dfdx(M*CellCount+i,((i-1)*(nelx*nely)+1):(i*(nelx*nely))) = 1000*(reshape(dv(1:nely,(((i-1)*nelx)+1):(i*nelx)),[1,nelx*nely])/(nelx*nely*volfrac_d(i)));
        end
        M = M + 1;
    end
    
    if LOCAL == 1
        for i = 1:CellCount
            fval(M*CellCount+i) = (sum(x_pde_hat_d{i}.^p))^(1/p)- vol_max_pNorm(i);
            dfdx(M*CellCount+i,((i-1)*(nelx*nely)+1):(i*(nelx*nely))) = TF'*(LF'\(LF\(TF*dfdx_pde_d{i}(:))));
            tmp = reshape(dfdx(M*CellCount+i,((i-1)*(nelx*nely)+1):(i*(nelx*nely))),[nely,nelx]);
            dfdx(M*CellCount+i,((i-1)*(nelx*nely)+1):(i*(nelx*nely))) = reshape(H_ind*(tmp(:).*dx_ind{i}(:)./Hs_ind),[1,nelx*nely]);
        end
        M = M + 1;
    end
    
    if ISO == 1
        for i = 1:CellCount
            fval(M*CellCount+i) = 10E6*(Q{i}(1,1)+Q{i}(2,2)-(Q{i}(1,2)+Q{i}(2,1))-4*Q{i}(3,3))^2 - epsilon_iso;
            dfdx(M*CellCount+i,((i-1)*(nelx*nely)+1):(i*(nelx*nely))) = 10E6*reshape((2*(Q{i}(1,1)+Q{i}(2,2)-(Q{i}(1,2)+Q{i}(2,1))-4*Q{i}(3,3))*(dQ{i}{1,1}+dQ{i}{2,2}-(dQ{i}{1,2}+dQ{i}{2,1})-4*dQ{i}{3,3})),[1,nelx*nely]);
        end
    end
    
    iter = loop;
    xval = reshape(x_cat,[CellCount*nelx*nely,1]);
    xmin = max(0.0,xval-move);
    xmax = min(1,xval+move);
    
    dfdx2 = 0*dfdx;
    a0 = 1;
    a = zeros(m,1);
    c_ = ones(m,1)*1000;
    d = zeros(m,1);
    mdof = 1:m;
    
    [xmma,ymma,zmma,lam,xsi,eta_,mu,zet,s,low,upp] = ...
        mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2,...
        f0val,df0dx,df0dx2,fval(1:m),dfdx(mdof,:),dfdx2(mdof,:),low,upp,a0,a,c_,d);
    
    xnew_cat = reshape(xmma,[nely,CellCount*nelx]);
    xold2 = xold1;
    xold1 = xval;
    
    if FilterDom == 1
        xTilde_cat(:) = (H*xnew_cat(:))./Hs;
    else
        for i = 1:CellCount
            xTilde_temp = xTilde_cat(1:nely,(((i-1)*nelx)+1):(i*nelx));
            xnew_temp = xnew_cat(1:nely,(((i-1)*nelx)+1):(i*nelx));
            xTilde_temp(:) = (H_ind*xnew_temp(:))./Hs_ind;
            xTilde_cat(1:nely,(((i-1)*nelx)+1):(i*nelx)) = xTilde_temp;
        end
    end
    
    xPhys_e_cat = (tanh(beta*eta) + tanh(beta*(xTilde_cat-eta))) / (tanh(beta*eta) + tanh(beta*(1-eta)));
    xPhys_i_cat = (tanh(beta*0.5) + tanh(beta*(xTilde_cat-0.5))) / (tanh(beta*0.5) + tanh(beta*(1-0.5)));
    xPhys_d_cat = (tanh(beta*eta) + tanh(beta*(xTilde_cat-eta))) / (tanh(beta*eta) + tanh(beta*(1-eta)));
    change = max(abs(xnew_cat(:)-x_cat(:)));
    x_cat = xnew_cat;
    
    for i = 1:CellCount
        x{i} = x_cat(1:nely,(((i-1)*nelx)+1):(i*nelx));
        xTilde{i} = xTilde_cat(1:nely,(((i-1)*nelx)+1):(i*nelx));
        xPhys_e{i} = xPhys_e_cat(1:nely,(((i-1)*nelx)+1):(i*nelx));
        xPhys_i{i} = xPhys_i_cat(1:nely,(((i-1)*nelx)+1):(i*nelx));
        xPhys_d{i} = xPhys_d_cat(1:nely,(((i-1)*nelx)+1):(i*nelx));
        
        x{i} = ApplySymmetry(x{i},sym);
        xTilde{i} = ApplySymmetry(xTilde{i},sym);
        xPhys_e{i} = ApplySymmetry(xPhys_e{i},sym);
        xPhys_i{i} = ApplySymmetry(xPhys_i{i},sym);
        xPhys_d{i} = ApplySymmetry(xPhys_d{i},sym);
    end
    
    i = CellCount + 1;
    l = 1;
    while i <= CaseCount
        j = 1;
        k = j + l;
        while k <= CellCount
            x{i} = horzcat(x{j},x{k});
            xTilde{i} = horzcat(xTilde{j},xTilde{k});
            xPhys_e{i} = horzcat(xPhys_e{j},xPhys_e{k});
            xPhys_i{i} = horzcat(xPhys_i{j},xPhys_i{k});
            xPhys_d{i} = horzcat(xPhys_d{j},xPhys_d{k});
            i  = i + 1;
            j = j + 1;
            k = j + l;
        end
        l = l + 1;
    end
    
    %% PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, ...
        mean(xPhys_i_cat(:)),change);
    
    %% UPDATE HEAVISIDE REGULARIZATION PARAMETER
    if beta < 256 && (loopbeta >= BetaUpdFreq || change <= 0.001)
        beta = BetaUpdFactor * beta;
        loopbeta = 0;
        change = 1;
        fprintf('Parameter beta increased to %g.\n',beta);
    end
    
    %% UPDATE ISOTROPY TOLERANCE PARAMETER
    if epsilon_iso > 0.001 && mod(loop,EpsilonUpdFreq) == 0
        epsilon_iso = EpsilonUpdFactor * epsilon_iso;
    end
    
    %% PLOT DENSITIES
    figure(1);
    colormap(gray); imagesc(1-xPhys_i_cat); caxis([0 1]); axis equal; axis off; drawnow;
    xPhys_cat = xPhys_i_cat;
    
    %% PLOT OBJECTIVE FUNCTION
    figure(2);
    clf;
    hold on
    for i = 1:CaseCount
        if i <= CellCount
            plot(c_hist(1:loop,i),'b')
        else
            plot(c_hist(1:loop,i),'g')
        end
    end
    % autoArrangeFigures(2,2,1);
end
end

%% FUNCTION: PREPARE FINITE ELEMENT ANALYSIS
function [nodenrs,edofMat,iK,jK] = PrepareFEA(nelx,nely)
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
end

%% FUNCTION: APPLY PERIODIC BOUNDARY CONDITIONS
function [U,ufixed,d1,d2,d3,d4,wfixed] = ApplyPeriodicBC(nelx,nely,nodenrs)
e0 = eye(3);
ufixed = zeros(8,3);
U = zeros(2*(nely+1)*(nelx+1),3);
alldofs = (1:2*(nely+1)*(nelx+1));
n1 = [nodenrs(end,[1,end]),nodenrs(1,[end,1])];
d1 = reshape([(2*n1-1);2*n1],1,8);
n3 = [nodenrs(2:end-1,1)',nodenrs(end,2:end-1)];
d3 = reshape([(2*n3-1);2*n3],1,2*(nelx+nely-2));
n4 = [nodenrs(2:end-1,end)',nodenrs(1,2:end-1)];
d4 = reshape([(2*n4-1);2*n4],1,2*(nelx+nely-2));
d2 = setdiff(alldofs,[d1 d3 d4]);
for j = 1:3
    ufixed(3:4,j) = [e0(1,j) e0(3,j)/2; e0(3,j)/2 e0(2,j)] * [nelx;0];
    ufixed(7:8,j) = [e0(1,j) e0(3,j)/2; e0(3,j)/2 e0(2,j)] * [0;nely];
    ufixed(5:6,j) = ufixed(3:4,j) + ufixed(7:8,j);
end
wfixed = [repmat(ufixed(3:4,:),nely-1,1); repmat(ufixed(7:8,:),nelx-1,1)];
end

%% FUNCTION: INITIALIZE DESIGN VARIABLES
function [x, xTilde, xPhys_e, xPhys_i, xPhys_d, xold1, xold2] = InitializeIteration(nelx,nely,beta,eta,volfrac,IC,x_imp,imp)
if imp == 0
    x = repmat(volfrac,nely,nelx);
    for i = 1:nelx
        for j = 1:nely
            if IC == 1
                if sqrt((i-nelx/2-0.5)^2+(j-nely/2-0.5)^2) < min(nelx,nely)/6
                    x(j,i) = volfrac/2;
                end
            elseif IC == 2
                if sqrt((i-nelx/4-0.5)^2+(j-nely/4-0.5)^2) < min(nelx,nely)/6
                    x(j,i) = volfrac/2;
                end
            elseif IC == 3
                if sqrt((i-nelx/4-0.5)^3+(j-nely/8-0.5)^2) < min(nelx/2,nely)/6
                    x(j,i) = volfrac/2;
                end
            end
        end
    end
else
    x = x_imp;
end

xTilde = x;
xPhys_e = (tanh(beta*eta) + tanh(beta*(xTilde-eta))) / (tanh(beta*eta) + tanh(beta*(1-eta)));
xPhys_i = (tanh(beta*0.5) + tanh(beta*(xTilde-0.5))) / (tanh(beta*0.5) + tanh(beta*(1-0.5)));
xPhys_d = (tanh(beta*eta) + tanh(beta*(xTilde-eta))) / (tanh(beta*eta) + tanh(beta*(1-eta)));

xold1 = reshape(x,[nely*nelx,1]);
xold2 = reshape(x,[nely*nelx,1]);
end

%% FUNCTION: FINITE ELEMENT ANALYSIS
function U = FEAnalysis(nelx,nely,KE,Emin,xPhys,penal,E0,iK,jK,d1,d2,d3,d4,ufixed,wfixed)
sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
K = sparse(iK,jK,sK); K = (K+K')/2;
Kr = [K(d2,d2), K(d2,d3)+K(d2,d4); K(d3,d2)+K(d4,d2), K(d3,d3)+K(d4,d3)+K(d3,d4)+K(d4,d4)];
U(d1,:) = ufixed;
U([d2,d3],:) = Kr\(-[K(d2,d1); K(d3,d1)+K(d4,d1)]*ufixed-[K(d2,d4); ...
    K(d3,d4)+K(d4,d4)]*wfixed);
U(d4,:) = U(d3,:) + wfixed;
end

%% FUNCTION: HOMOGENIZATION
function [Q, dQ] = Homogenization(nelx,nely,U,edofMat,Emin,E0,xPhys,penal,KE)
qe = cell(3,3);
Q = zeros(3);
dQ = cell(3,3);
for i = 1:3
    for j = 1:3
        U1 = U(:,i); U2 = U(:,j);
        qe{i,j} = reshape(sum((U1(edofMat)*KE).*U2(edofMat),2),...
            nely,nelx)/(nelx*nely);
        Q(i,j) = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*qe{i,j}));
        dQ{i,j} = penal*(E0-Emin)*xPhys.^(penal-1).*qe{i,j};
    end
end
end

%% FUNCTION: OBJECTIVE
function [c,dc] = Objective(Q,dQ,OBJECTIVE,loop)
if OBJECTIVE == 1
    c = -(Q(1,1)+Q(2,2)+Q(1,2)+Q(2,1));
    dc = -(dQ{1,1}+dQ{2,2}+dQ{1,2}+dQ{2,1});
elseif OBJECTIVE == 2
    c = -Q(3,3);
    dc = -dQ{3,3};
elseif OBJECTIVE == 3
    c = Q(1,2) - (0.8^loop)*(Q(1,1)+Q(2,2));
    dc = dQ{1,2} -(0.8^loop)*(dQ{1,1}+dQ{2,2});
end
end

%% APPLY SYMMETRY
function x = ApplySymmetry(x,sym)
nelx = size(x,2);

% 4-fold symmetry across diagonals
if sym == 1
    x = (x+x')/2;
    x = rot90(x);
    x = (x+x')/2;
    x = rot90(x);
    
    % 4-fold symmetry across normals
elseif sym == 2
    for row = 1:nelx
        for col = 1:nelx/2
            x(row,col) = (x(row,col) + x(row,nelx-col+1))/2;
            x(row,nelx+1-col) = x(row,col);
        end
    end
    
    for col = 1:nelx
        for row = 1:nelx/2
            x(row,col) = (x(row,col) + x(nelx-row+1,col))/2;
            x(nelx+1-row,col) = x(row,col);
        end
    end
    
    % 4-fold rotational symmetry
elseif sym == 3
    A = x;
    A90 = rot90(A);
    A180 = rot90(A,2);
    A270 = rot90(A,3);
    E = zeros(nelx);
    temp = 0.25*(A(1:nelx/2,1:nelx/2)+A90(1:nelx/2,1:nelx/2)+A180(1:nelx/2,1:nelx/2)+A270(1:nelx/2,1:nelx/2));
    for ii = 1:3
        E(1:nelx/2,1:nelx/2) = temp;
        E = rot90(E);
    end
    E(1:nelx/2,1:nelx/2) = temp;
    x = E;
    
    % Cubic symmetry
elseif sym == 4
    A = x;
    A90 = rot90(A);
    A180 = rot90(A,2);
    A270 = rot90(A,3);
    E = zeros(nelx);
    temp = 0.25*(A(1:nelx/2,1:nelx/2)+A90(1:nelx/2,1:nelx/2)+A180(1:nelx/2,1:nelx/2)+A270(1:nelx/2,1:nelx/2));
    for ii = 1:3
        E(1:nelx/2,1:nelx/2) = temp;
        E = rot90(E);
    end
    E(1:nelx/2,1:nelx/2) = temp;
    x = E;
    
    for row = 1:nelx
        for col = 1:nelx/2
            x(row,col) = (x(row,col) + x(row,nelx-col+1))/2;
            x(row,nelx+1-col) = x(row,col);
        end
    end
    
    for col = 1:nelx
        for row = 1:nelx/2
            x(row,col) = (x(row,col) + x(nelx-row+1,col))/2;
            x(nelx+1-row,col) = x(row,col);
        end
    end
end
end

% Publication
% Eric Garner, Helena Kolken, Charlie C.L. Wang, Amir A. Zadpoor, Jun Wu
% Compatibility in Microstructural Optimization for Additive Manufacturing
% Additive Manufacturing, submitted 2018

% The code was developed based on the periodic material microstructure design
% code, by L. Xia and P. Breitkopf, 2014
