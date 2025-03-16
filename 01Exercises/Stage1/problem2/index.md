<script type="text/javascript" async
  src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js">
</script>
# Problem 2: Implement other boundary conditions

<p align="center">
  <figure align="center">
    <img src="../../image/stage1/problem2/Problem2.png" width="90%">
  </figure>
</p>

- Design Domain: A 2D rectangular plane, discretized into a grid of nelx=120 by nely=20 elements.
- Design Variables, objective function,  constraint are the same as Problem 1.

- Optimization parameters: The volume fraction is set to 0.5, the penalization factor is set to 3, the filter radius is set to 4.8, and the filter type is set to sensitivity filter.
## Two simultaneous point loads
- Loading Conditions: Two downward point loads are applied at the top edge of the design domain, specifically at the one-quarter and three-quarters positions. We compared the cases where the forces are applied upwards versus downwards.
```matlab
F(2*(nelx/4)*(nely+1)+2,1)=1; F(2*(3*nelx/4)*(nely+1)+2,1)=1;
```
### Case1
- Boundary Conditions: The nodes at the bottom left and right corners of the design domain are fully fixed, meaning that both the x and y degrees of freedom at these nodes are constrained.
```matlab
fixeddofs = [2*(nely+1)-1, 2*(nely+1), 2*(nelx+1)*(nely+1)-1, 2*(nelx+1)*(nely+1)];
```
<p align="center">
  <figure align="center">
    <img src="../../image/stage1/problem2/result_two_point_loads_both_xy_iter8_obj82.6921.png" width="90%">
  </figure>
</p>



### Case2
- Boundary Conditions: The nodes at the bottom left and right corners of the design domain are partially fixed, meaning that only the y degrees of freedom at these nodes are constrained, while the x degrees of freedom remain free.
```matlab
fixeddofs = [2*(nely+1), 2*(nelx+1)*(nely+1)];
```
<p align="center">
  <figure align="center">
    <img src="../../image/stage1/problem2/result_two_point_loads_only_y_iter33_obj226.4264.png" width="90%">
  </figure>
</p>



## Distributed load
- Loading Conditions: A series of vertical forces are applied along the top edge of the design domain.
```matlab
for i = 1:nelx+1
    F(2*(i-1)*(nely+1)+2, 1) = -1;
end
```
### Case1
- Boundary Conditions: The nodes at the bottom left and right corners of the design domain are fully fixed, meaning that both the x and y degrees of freedom at these nodes are constrained.
```matlab
fixeddofs = [2*(nely+1)-1, 2*(nely+1), 2*(nelx+1)*(nely+1)-1, 2*(nelx+1)*(nely+1)];
```
<p align="center">
  <figure align="center">
    <img src="../../image/stage1/problem2/result_distributed_load_both_xy_iter13_obj277578.1956.png" width="90%">
  </figure>
</p>


### Case2
- Boundary Conditions: The nodes at the bottom left and right corners of the design domain are partially fixed, meaning that only the y degrees of freedom at these nodes are constrained, while the x degrees of freedom remain free.
```matlab
fixeddofs = [2*(nely+1), 2*(nelx+1)*(nely+1)];
```
<p align="center">
  <figure align="center">
    <img src="../../image/stage1/problem2/result_distributed_load_only_y_iter18_obj867218.2148.png" width="90%">
  </figure>
</p>


## Comments and Conclusion
- Experiments show that the direction of the force does not change the optimization results (identical result images are not displayed in the report to save space). In fact, from the perspective of theoretical derivation, \\(KU = F\\) is a linear equation. When \\(F\\) changes to \\(-F\\), \\(U\\) changes to \\(-U\\). However, the objective function \\(c\\) is a quadratic function of \\(U\\), so changing \\(U\\) to \\(-U\\) does not affect the value of \\(c\\), and thus does not influence the update process of \\(x\\) at each step.
- We applied a downward force component to each node in the first row under the distributed load condition. **Should we expect the material density at the points where the force is applied to be 1 in the optimization results?**
- In this problem, I chose the filter radius as rmin = nelx * 0.04 = 120 * 0.04 = 4.8 based on the experience from Assignment 1. However, I overlooked the fact that in this problem, I did not reduce the problem to half of the material based on symmetry, so the ratio of nelx to nely is 6, not 3. As a result, the filter radius rmin = 4.8 might be too large for nely, which could have led to large gray areas in the optimization results.
- In this problem, I discussed two boundary conditions: one where only the y components at the bottom-left and bottom-right corners are fixed, and another where both the x and y components are fixed. I'm not sure which boundary condition better reflects the actual physics. **How should we determine the boundary conditions in real-world problems?**

## Source code
```matlab
function problem2()
    %% FIXED PARAMETERS
    nelx = 120;      % Fixed number of elements along x-axis
    nely = 20;       % Fixed number of elements along y-axis
    volfrac = 0.5;   % Fixed volume fraction
    penal = 3;       % Fixed penalization factor
    rmin = 4.8;      % Fixed filter radius
    ft = 1;          % Fixed filter type (sensitivity filter)

    %% MATERIAL PROPERTIES
    E0 = 1;
    Emin = 1e-9;
    nu = 0.3;

    %% LOAD AND BOUNDARY CONDITIONS
    load_conditions = {'two_point_loads', 'distributed_load'};
    boundary_conditions = {'only_y', 'both_xy'};
    
    for lc = 1:2
        for bc = 1:2
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

            %% DEFINE LOADS
            F = sparse(2*(nely+1)*(nelx+1), 1);
            U = zeros(2*(nely+1)*(nelx+1),1);
            if strcmp(load_conditions{lc}, 'two_point_loads')
                F(2*(nelx/4)*(nely+1)+2,1) = -1;
                F(2*(3*nelx/4)*(nely+1)+2,1) = -1;
            elseif strcmp(load_conditions{lc}, 'distributed_load')
                for i = 1:nelx+1
                    F(2*(i-1)*(nely+1)+2, 1) = 1;
                end
            end

            %% DEFINE SUPPORTS
            if strcmp(boundary_conditions{bc}, 'only_y')
                fixeddofs = [2*(nely+1), 2*(nelx+1)*(nely+1)];
            elseif strcmp(boundary_conditions{bc}, 'both_xy')
                fixeddofs = [2*(nely+1)-1, 2*(nely+1), 2*(nelx+1)*(nely+1)-1, 2*(nelx+1)*(nely+1)];
            end

            alldofs = 1:2*(nely+1)*(nelx+1);
            freedofs = setdiff(alldofs, fixeddofs);

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
            loop = 0;
            change = 1;
            c_history = [];

            %% START ITERATION
            while change > 0.01
                loop = loop + 1;

                %% FE-ANALYSIS
                sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
                K = sparse(iK,jK,sK); K = (K+K')/2;
                U(freedofs) = K(freedofs,freedofs)\F(freedofs);

                %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
                ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
                c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
                dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
                dv = ones(nely,nelx);
                c_history = [c_history c];  % Store objective function value

                %% FILTERING/MODIFICATION OF SENSITIVITIES
                if ft == 1
                    dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
                elseif ft == 2
                    dc(:) = H*(dc(:)./Hs);
                    dv(:) = H*(dv(:)./Hs);
                end

                %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
                l1 = 0; l2 = 1e9; move = 0.2;
                while (l2-l1)/(l1+l2) > 1e-3
                    lmid = 0.5*(l2+l1);
                    xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dc./dv/lmid)))));
                    if ft == 1
                        xPhys = xnew;
                    elseif ft == 2
                        xPhys(:) = (H*xnew(:))./Hs;
                    end
                    if sum(xPhys(:)) > volfrac*nelx*nely, l1 = lmid; else l2 = lmid; end
                end
                change = max(abs(xnew(:)-x(:)));
                x = xnew;
            end

           %% PLOT DENSITIES AND OBJECTIVE FUNCTION HISTORY
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
save_filename = sprintf('result_%s_%s_iter%d_obj%.4f.png', load_conditions{lc}, boundary_conditions{bc}, loop, c);
saveas(gcf, save_filename);
close all;
        end
    end
end
```

