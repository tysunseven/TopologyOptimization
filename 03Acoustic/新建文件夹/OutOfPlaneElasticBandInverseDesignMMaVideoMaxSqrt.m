function [xfinal, figure, loop, time, bdmax, bdmin, bdmax_x, bdmax_y, bdmin_x, bdmin_y, bdgap,...
    gbmax, gbmin, gbmax_x, gbmax_y, gbmin_x, gbmin_y, gbgap]...
    =OutOfPlaneElasticBandInverseDesignMMaVideoMaxSqrt(x,mode,num_random_points)
Lx = 2; nelx = size(x,2); nely = nelx; h = Lx/nelx; num_modes = 5; p = 48; q = -p;
rho1 = 1; rho2 = 2; E1 = 4; E2 = 20; nu = 0.34; mu1 = E1/(2*(1+nu)); mu2 = E2/(2*(1+nu));
kpoint = struct( 'mu_x', 0, 'mu_y', 0 );
[boundary_mu_x, boundary_mu_y, path_distance,xtick_pos, xtick_labels] = generate_band_path_OABCOB();
kpoints = repmat(kpoint, numel(boundary_mu_x)+num_random_points,1);
for i = 1:numel(boundary_mu_x)
    kpoints(i).mu_x = boundary_mu_x(i); kpoints(i).mu_y = boundary_mu_y(i);
end
eigenvalues = zeros(length(boundary_mu_x)+num_random_points, num_modes);
[edofMat, Ke, Me, iIndex, jIndex] = init_fem(nelx,nely,h); [row, col, fixT] = init_trans(nelx,nely);
loop = 0; change = 1; [m, n, xold1, xold2, low, upp, a0, a, c_MMA, d, move] = init_mma(nelx, nely, x, 1e-3);
figHandles = initOptimizationFigures(x, num_modes, xtick_pos, xtick_labels);
max_vals = []; min_vals = []; gap_vals = []; approxmax_vals = []; approxmin_vals = []; approxgap_vals = [];
%% 循环体
tic;
while change>0.01 || loop<20
    loop = loop + 1; xmin = max( x(:) - move, 0 ); xmax = min( x(:) + move, 1 );

    mu = mu1*(1-x)+mu2*x;
    rho = rho1*(1-x)+rho2*x;
    sKa = reshape(Ke(:)*(mu(:))',16*nelx*nely,1);
    sMa = reshape(Me(:)*(rho(:))',16*nelx*nely,1);
    K = sparse(iIndex,jIndex,sKa); K = (K+K')/2;
    M = sparse(iIndex,jIndex,sMa); M = (M+M')/2;
    
    random_mu_x = pi * rand(num_random_points, 1); random_mu_y = pi * rand(num_random_points, 1);
    if num_random_points > 0
        for i = numel(boundary_mu_x)+1:numel(boundary_mu_x)+numel(num_random_points)
        kpoints(i).mu_x = random_mu_x(i-numel(boundary_mu_x)); kpoints(i).mu_y = random_mu_y(i-numel(boundary_mu_x));
        end
    end

    dcmax = zeros(nely, nelx); dcmin = zeros(nely, nelx);

    parfor i = 1:length(kpoints)
        T = create_T(kpoints(i).mu_x, kpoints(i).mu_y, nelx, nely, row, col, fixT);
        K_tilde = T' * K * T; M_tilde = T' * M * T;
        [V, D] = eigs(K_tilde, M_tilde, num_modes, 'sm');
        eigenvalues(i,:) = sort(sqrt(abs(real(diag(D)))));
        phi = T*V(:,mode);  psi = T*V(:,mode+1);
        ceKmax = (mu2-mu1)*reshape(sum((phi(edofMat)*Ke).*phi(edofMat),2),nely,nelx);
        ceMmax = (rho2-rho1)*reshape(sum((phi(edofMat)*Me).*phi(edofMat),2),nely,nelx);
        dki1dx = real((ceKmax-D(mode,mode)*ceMmax)/(2*sqrt(real(D(mode,mode)))*phi'*M*phi));
        dcmax = dcmax + sqrt(real(D(mode,mode)))^(p-1)*dki1dx;
        ceKmin = (mu2-mu1)*reshape(sum((psi(edofMat)*Ke).*psi(edofMat),2),nely,nelx);
        ceMmin = (rho2-rho1)*reshape(sum((psi(edofMat)*Me).*psi(edofMat),2),nely,nelx);
        dki2dx = real((ceKmin-D(mode+1,mode+1)*ceMmin)/(2*sqrt(real(D(mode+1,mode+1)))*psi'*M*psi));
        dcmin = dcmin + sqrt(real(D(mode+1,mode+1)))^(q-1)*dki2dx;
    end
    
    dcmax = dcmax * sum(eigenvalues(:, mode).^p)^(1/p-1);
    dcmin = dcmin * sum(eigenvalues(:, mode+1).^q).^(1/q-1);
    apxmax = sum(eigenvalues(:, mode).^p).^(1/p);
    apxmin = sum(eigenvalues(:, mode+1).^q).^(1/q);
    c = apxmax-apxmin+20;
    dc = dcmax-dcmin;
    % MMA优化
    
    f0val = c;                                    % scalar
    df0dx = dc(:);                                % column 
    dv = ones(nely,nelx); 
    fval  = (sum(x(:))/(1000*n) - 1);            % scalar volumn constrain
    dfdx  = dv(:)'/ (0.5*n);                  % row


    [ xmma, ~, ~, ~, ~, ~, ~, ~, ~, low, upp ] = ...
        mmasub( m, n, loop, x(:), xmin, xmax, xold1, xold2, f0val,df0dx,fval,dfdx,low,upp,a0,a,c_MMA,d);
    xold2 = xold1(:); xold1 = x(:); x = reshape(xmma,nely,nelx); change = max(abs(x(:)-xold1(:)));
      

    [figHandles, max_vals, min_vals, approxmax_vals, approxmin_vals, approxgap_vals, gap_vals] = ...
    updateFigures(figHandles, eigenvalues, x, mode, loop,...
                              path_distance, boundary_mu_x, num_modes,...
                              apxmax, apxmin, max_vals, min_vals,...
                              approxmax_vals, approxmin_vals, approxgap_vals, gap_vals);

    fprintf('Iter:%4i max:%7.4f apxmax:%7.4f min:%7.4f apxmin:%7.4f gap:%7.4f apxgap:%7.4f chg:%5.3f\n', ...
        loop, max(eigenvalues(:,mode)), apxmax, min(eigenvalues(:,mode+1)), apxmin,  min(eigenvalues(:,mode+1))-max(eigenvalues(:,mode)),apxmin-apxmax,  change);
end
time = toc;
xfinal = x;
figure = figHandles.fig;
[bdmax, bdmin, bdmax_x, bdmax_y, bdmin_x, bdmin_y, bdgap,...
    gbmax, gbmin, gbmax_x, gbmax_y, gbmin_x, gbmin_y, gbgap]...
    =postFigures(mode,num_modes,nelx,nely,row,col,fixT,K,M);