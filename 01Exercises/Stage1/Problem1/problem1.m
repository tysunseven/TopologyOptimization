function problem1()
    % 初始化Excel数据存储
    results = {};
    results{1, 1} = 'nelx';
    results{1, 2} = 'nely';
    results{1, 3} = 'volfrac';
    results{1, 4} = 'penal';
    results{1, 5} = 'rmin';
    results{1, 6} = 'ft';
    results{1, 7} = 'loop';
    results{1, 8} = 'c';
    
    % 记录结果的行数索引
    row_index = 2;

    % 第一层循环遍历nelx的值
    % for nelx = [60, 120, 240]
    for nelx = [60]
        nely = nelx / 3;
        
        % 第二层循环遍历volfrac的值
        for volfrac = [0.5]
            
            % 第三层循环遍历penal的值
            for penal = [3]
                
                % 第四层循环遍历rmin的值
                for rmin = [2.4]
                    
                    % 第五层循环遍历ft的值
                    for ft = [0, 1]
                        
                        % 调用top88函数并返回loop和c
                        [loop, c] = top88(nelx, nely, volfrac, penal, rmin, ft);
                        
                        % 将参数和结果存储在results中
                        results{row_index, 1} = nelx;
                        results{row_index, 2} = nely;
                        results{row_index, 3} = volfrac;
                        results{row_index, 4} = penal;
                        results{row_index, 5} = rmin;
                        results{row_index, 6} = ft;
                        results{row_index, 7} = loop;
                        results{row_index, 8} = c;
                        
                        % 更新行索引
                        row_index = row_index + 1;
                    end
                end
            end
        end
    end
    
    % 将结果保存到Excel文件中
    result_table = cell2table(results(2:end,:), 'VariableNames', results(1,:));
    writetable(result_table, 'Problem1_Results.xlsx');
end

function [loop, c] = top88(nelx,nely,volfrac,penal,rmin,ft)
%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;
nu = 0.3;
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
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F = sparse(2,1,-1,2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);
fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
alldofs = [1:2*(nely+1)*(nelx+1)];
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
loop = 0;
change = 1;
objective_history = []; % Initialize the array to store objective function values
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
  objective_history = [objective_history, c]; % Store current objective function value
  %% FILTERING/MODIFICATION OF SENSITIVITIES
  if ft==0
      % No filtering if ft == 0
  elseif ft == 1
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
    if ft == 0
      xPhys = xnew;
    elseif ft == 1
      xPhys = xnew;
    elseif ft == 2
      xPhys(:) = (H*xnew(:))./Hs;
    end
    if sum(xPhys(:)) > volfrac*nelx*nely, l1 = lmid; else l2 = lmid; end
  end
  change = max(abs(xnew(:)-x(:)));
  x = xnew;
  %% CHECK MAXIMUM NUMBER OF ITERATIONS
  if loop >= 2000
      fprintf('Max iterations reached.\n');
      break;
  end
  %% PRINT RESULTS
  fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c, ...
    mean(xPhys(:)),change);
  %% PLOT DENSITIES
  colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
end
%% SAVE FINAL STRUCTURE IMAGE
structure_filename = sprintf('Problem1structure_nelx%d_nely%d_volfrac%.2f_penal%.2f_rmin%.2f_ft%d_loop%d_c%.4f.png', ...
    nelx, nely, volfrac, penal, rmin, ft, loop, c);
saveas(gcf, structure_filename);

%% PLOT AND SAVE OBJECTIVE FUNCTION HISTORY
figure;
plot(1:loop, objective_history, '-o', 'LineWidth', 2);
xlabel('Iteration');
ylabel('Objective Function Value');
title('Objective Function vs Iteration');
grid on;
objective_filename = sprintf('Problem1objective_nelx%d_nely%d_volfrac%.2f_penal%.2f_rmin%.2f_ft%d_loop%d_c%.4f.png', ...
    nelx, nely, volfrac, penal, rmin, ft, loop, c);
saveas(gcf, objective_filename);
end