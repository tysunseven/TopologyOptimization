function OutOfPlaneElasticBandInverseDesignMMa(x)
nelx=size(x,1);
nely=size(x,2);
Lx=2;
h=Lx/nelx;
[space1,space2,time1,time2] = outofplaneconst();

num_modes = 5; num_random_points = 0;

kpoint = struct( 'mu_x', 0, 'mu_y', 0 );

[boundary_mu_x, boundary_mu_y, path_distance] = generate_band_path_OABCOB();

kpoints = repmat(kpoint, numel(boundary_mu_x)+num_random_points,1);
for i = 1:numel(boundary_mu_x)
    kpoints(i).mu_x = boundary_mu_x(i); kpoints(i).mu_y = boundary_mu_y(i);
end

eigenvalues = zeros(length(boundary_mu_x)+num_random_points, num_modes);
eigenvectors = zeros(length(boundary_mu_x)+num_random_points, (nelx+1)*(nely+1), num_modes);


[edofMat, Ke, Me, iIndex, jIndex] = init_fem(nelx,nely,h); [row, col, fixT] = init_trans(nelx,nely);
[m, n, xmin, xmax, xold1, xold2, low, upp, a0, a, c_MMA, d, move] = init_mma(nelx, nely, x, 1000);
[fig, ax_band, ax_density, band_plots, density_plot] = init_plots(x, num_modes);

loop = 0; change = 1;
while change>0.01
    loop = loop + 1; xmin = max( x(:) - move, 0 ); xmax = min( x(:) + move, 1 );
    [K, M] = loop_fem(space1,space2,time1,time2,x,iIndex,jIndex,Ke,Me,nelx,nely);

    random_mu_x = pi * rand(num_random_points, 1); random_mu_y = pi * rand(num_random_points, 1);
    if num_random_points > 0
        for i = numel(boundary_mu_x)+1:numel(boundary_mu_x)+numel(num_random_points)
        kpoints(i).mu_x = random_mu_x(i-numel(boundary_mu_x)); kpoints(i).mu_y = random_mu_y(i-numel(boundary_mu_x));
        end
    end

    parfor i = 1:length(kpoints)
        T = create_T(kpoints(i).mu_x, kpoints(i).mu_y, nelx, nely, row, col, fixT);
        K_tilde = T' * K * T; M_tilde = T' * M * T;
        [V, D] = eigs(K_tilde, M_tilde, num_modes, 'sm');
        eigenvalues(i,:) = sort(sqrt(abs(real(diag(D)))));
        eigenvectors(i,:,:) = T * V;
    end
    
    [f0val, df0dx, fval, dfdx] = min_single(...
        40,kpoints, eigenvalues, eigenvectors, space2, space1, time2, time1, edofMat, Ke, Me, M, x, n);

    [ xmma, ~, ~, ~, ~, ~, ~, ~, ~, low, upp ] = ...
        mmasub( m, n, loop, x(:), xmin, xmax, xold1, xold2, f0val,df0dx,fval,dfdx,low,upp,a0,a,c_MMA,d);
    xold2 = xold1(:); xold1 = x(:); x = reshape(xmma,nely,nelx); change = max(abs(x(:)-xold1(:)));
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,f0val,mean(x(:)),change);
    

    for l = 1:num_modes
        set(band_plots(l), 'XData', path_distance, 'YData', eigenvalues(1:numel(boundary_mu_x), l));
    end
    
    set(density_plot, 'CData', real(1 - x)); drawnow limitrate;
end