% Efficient code to do optimization using rank-2 microstructures, Feb. 2021
% variables:
% mu: Relative layer widths
% a: layer orietnation of first layer
%
% Coordinate system is defined as y positive downwards.
% Angles are flipped for plotting due to imagesc coordinate system

function topRank2(nelX,nelY,rMin,volFrac)
%% Material properties
E = 1;                  % Young's modulus
nu = 0.3;           % Poisson's ratio
EMin = 1e-1;    % Young's modulus of isotropic background material (for stability)
EMinFinal = 1e-6; % Young's modulus of isotropic background material (for stability)
%% Prepare finite element analysis
KE011 = [1/3,0,-1/3,0,-1/6,0,1/6,0;zeros(1,8); -1/3,0,1/3,0,1/6,0,-1/6,0;zeros(1,8);-1/6,0,1/6,0,1/3,0,-1/3,0;zeros(1,8); 1/6,0,-1/6,0,-1/3,0,1/3,0;zeros(1,8)];
KE012 = [0,1/4,0,1/4,0,-1/4,0,-1/4;1/4,0,-1/4,0,-1/4,0,1/4,0;0,-1/4,0,-1/4,0,1/4,0,1/4;1/4,0,-1/4,0,-1/4,0,1/4,0; 0,-1/4,0,-1/4,0,1/4,0,1/4;-1/4,0,1/4,0,1/4,0,-1/4,0;0,1/4,0,1/4,0,-1/4,0,-1/4;-1/4,0,1/4,0,1/4,0,-1/4,0];
KE013 = [1/2,1/3,0,-1/3,-1/2,-1/6,0,1/6;1/3,0,-1/3,0,-1/6,0,1/6,0;0,-1/3,-1/2,1/3,0,1/6,1/2,-1/6;-1/3,0,1/3,0,1/6,0,-1/6,0;-1/2,-1/6,0,1/6,1/2,1/3,0,-1/3;-1/6,0,1/6,0,1/3,0,-1/3,0;0,1/6,1/2,-1/6,0,-1/3,-1/2,1/3;1/6,0,-1/6,0,-1/3,0,1/3,0];
KE022 = [zeros(1,8);0,1/3,0,1/6,0,-1/6,0,-1/3;zeros(1,8);0,1/6,0,1/3,0,-1/3,0,-1/6;zeros(1,8);0,-1/6,0,-1/3,0,1/3,0,1/6;zeros(1,8);0,-1/3,0,-1/6,0,1/6,0,1/3];
KE023 = [0,1/3,0,1/6,0,-1/6,0,-1/3;1/3,1/2,1/6,0,-1/6,-1/2,-1/3,0;0,1/6,0,1/3,0,-1/3,0,-1/6;1/6,0,1/3,-1/2,-1/3,0,-1/6,1/2;0,-1/6,0,-1/3,0,1/3,0,1/6;-1/6,-1/2,-1/3,0,1/3,1/2,1/6,0;0,-1/3,0,-1/6,0,1/6,0,1/3;-1/3,0,-1/6,1/2,1/6,0,1/3,-1/2];
KE033 = [1/3,1/4,1/6,-1/4,-1/6,-1/4,-1/3,1/4;1/4,1/3,1/4,-1/3,-1/4,-1/6,-1/4,1/6;1/6,1/4,1/3,-1/4,-1/3,-1/4,-1/6,1/4;-1/4,-1/3,-1/4,1/3,1/4,1/6,1/4,-1/6;-1/6,-1/4,-1/3,1/4,1/3,1/4,1/6,-1/4;-1/4,-1/6,-1/4,1/6,1/4,1/3,1/4,-1/3;-1/3,-1/4,-1/6,1/4,1/6,1/4,1/3,-1/4;1/4,1/6,1/4,-1/6,-1/4,-1/3,-1/4,1/3];
BMatrixMid = 1/8*[-1,0,1,0,1,0,-1,0;0,-1,0,-1,0,1,0,1;-1,-1,-1,1,1,1,1,-1];
nodeNumbers = reshape(1:(1+nelX)*(1+nelY),1+nelY,1+nelX);
edofVec = reshape(2*nodeNumbers(1:end-1,1:end-1)+1,nelX*nelY,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nelY+[2 3 0 1] -2 -1],nelX*nelY,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelX*nelY,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelX*nelY,1);
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F = sparse(2,1,-1,2*(nelY+1)*(nelX+1),1);
U = zeros(2*(nelY+1)*(nelX+1),1);
stress = zeros(nelX*nelY,3);
fixedDofs = union([1:2:2*(nelY+1)],[2*(nelX+1)*(nelY+1)]);
allDofs = [1:2*(nelY+1)*(nelX+1)];
freeDofs = setdiff(allDofs,fixedDofs);
%% PREPARE FILTER
iH = ones(nelX*nelY*(2*(ceil(rMin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelX
    for j1 = 1:nelY
        element1 = (i1-1)*nelY+j1;
        for i2 = max(i1-(ceil(rMin)-1),1):min(i1+(ceil(rMin)-1),nelX)
            for j2 = max(j1-(ceil(rMin)-1),1):min(j1+(ceil(rMin)-1),nelY)
                element2 = (i2-1)*nelY+j2;
                k = k+1;
                iH(k) = element1;
                jH(k) = element2;
                sH(k) = max(0,rMin-sqrt((i1-i2)^2+(j1-j2)^2));
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% INITIALIZE ITERATION
mu = repmat(1-(sqrt(1-volFrac)),nelX*nelY,2);
muPhys = mu; muNew = mu;
a = zeros(nelX*nelY,1);
loop = 0;
change = [1,1];
stopCrit = 0.01;
loopMax = 400;
%% Start iterations
while max(change) > stopCrit && loop < loopMax
    loop = loop + 1;
    %% FE-analysis
    [EH11,EH12,EH13,EH22,EH23,EH33,~,~,~,~,~,~,~,~,~,~,~,~] = getElasticityDerivativesRank2(a,muPhys,E,nu,EMin);
    sK = KE011(:)*EH11'+KE012(:)*EH12'+KE013(:)*EH13'+KE022(:)*EH22'+KE023(:)*EH23'+KE033(:)*EH33';
    K = sparse(iK,jK,sK(:));
    K = (K+K')/2;
    U(freeDofs) = K(freeDofs,freeDofs)\F(freeDofs);
    %% ANGLE UPDATE
    Ue = U(edofMat');
    strain = BMatrixMid*Ue;
    stress(:,1) = EH11.*strain(1,:)'+EH12.*strain(2,:)'+EH13.*strain(3,:)';
    stress(:,2) = EH12.*strain(1,:)'+EH22.*strain(2,:)'+EH23.*strain(3,:)';
    stress(:,3) = EH13.*strain(1,:)'+EH23.*strain(2,:)'+EH33.*strain(3,:)';
    aNew(:) = atan2(2*stress(:,3),(stress(:,1)-stress(:,2)))/2;
    change(2) = max(abs(a(:)-aNew(:)));
    a(:) = aNew(:);
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    [EH11,EH12,EH13,EH22,EH23,EH33,dEH11dMu1,dEH12dMu1,dEH13dMu1,dEH22dMu1,dEH23dMu1,dEH33dMu1,dEH11dMu2,dEH12dMu2,dEH13dMu2,dEH22dMu2,dEH23dMu2,dEH33dMu2] = getElasticityDerivativesRank2(a,muPhys,E,nu,EMin);
    Ue2 =  reshape(bsxfun(@times,reshape(Ue,[1 8 nelX*nelY]),reshape(Ue,[8 1 nelX*nelY])),64,nelX*nelY);
    sK = KE011(:)*EH11'+KE012(:)*EH12'+KE013(:)*EH13'+KE022(:)*EH22'+KE023(:)*EH23'+KE033(:)*EH33';
    dsKdMu1 = KE011(:)*dEH11dMu1'+KE012(:)*dEH12dMu1'+KE013(:)*dEH13dMu1'+KE022(:)*dEH22dMu1'+KE023(:)*dEH23dMu1'+KE033(:)*dEH33dMu1';
    dsKdMu2 = KE011(:)*dEH11dMu2'+KE012(:)*dEH12dMu2'+KE013(:)*dEH13dMu2'+KE022(:)*dEH22dMu2'+KE023(:)*dEH23dMu2'+KE033(:)*dEH33dMu2';
    ce = reshape(sum(sK.*Ue2,1),nelX*nelY,1);
    c = sum(ce,1);
    dcdMu1 = -reshape(sum(dsKdMu1.*Ue2,1),nelX*nelY,1);
    dcdMu2 = -reshape(sum(dsKdMu2.*Ue2,1),nelX*nelY,1);
    dvdMu1 = (1-muPhys(:,2));
    dvdMu2 = (1-muPhys(:,1));
    %% FILTERING/MODIFICATION OF SENSITIVITIES
    dcdMu1(:) = H*(dcdMu1(:)./Hs);
    dcdMu2(:) = H*(dcdMu2(:)./Hs);
    dvdMu1(:) = H*(dvdMu1(:)./Hs);
    dvdMu2(:) = H*(dvdMu2(:)./Hs);
    %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
    l1 = 0; l2 = 1e5; moveMu = 0.1;
    while (l2-l1)/(l1+l2) > 1e-3
        lmid = 0.5*(l2+l1);
        muNew = max(0,max(mu-moveMu,min(1.,min(mu+moveMu,(mu.*sqrt(-([dcdMu1,dcdMu2])./lmid./([dvdMu1,dvdMu2])))))));
        muPhys(:,1) = (H*muNew(:,1))./Hs;
        muPhys(:,2) = (H*muNew(:,2))./Hs;
        rho = muPhys(:,1)+muPhys(:,2)-muPhys(:,1).*muPhys(:,2);
        if sum(rho,1) >volFrac*nelX*nelY, l1 = lmid; else l2 = lmid; end
    end
    change(1) = max(abs(muNew(:)-mu(:)));
    mu(:) = muNew(:);
    %% Continuation on background stiffness
    if  (EMin > (1.001*EMinFinal) ) && (mod(loop,50) == 0 || max(change) < stopCrit)
        EMin = 0.1*EMin; change = [1,1];
        fprintf(' Reduced background stiffness to .:%8.7f\n',EMin);
    end
    %% PRINT RESULTS and plot
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f ch1.:%7.3f ch2.:%7.3f  \n',loop,c, ...
        sum(rho,1)/(nelX*nelY),change(1),change(2));
    %% PLOT DENSITIES
    figure(1);clf;colormap(gray); imagesc(1-reshape(rho,nelY, nelX)); caxis([0 1]); axis equal; axis off; drawnow;
end
%% Plot final design
figure(2);clf;
plotFinalDesign(nelX,nelY,rho,muPhys,a)
end

function [EH11,EH12,EH13,EH22,EH23,EH33] = rotateLocaToGlobal(a,EL11,EL12,EL22,EL33)
EH11 = (cos(a).^4).*(EL11-2*EL12+EL22-4*EL33)+(cos(a).^2).*(2*EL12-2*EL22+4*EL33)+EL22;
EH12 = -(cos(a).^4).*(EL11-2*EL12+EL22-4*EL33)+(cos(a).^2).*(EL11-2*EL12+EL22-4*EL33)+EL12;
EH13 = sin(a).*cos(a).*((EL11-2*EL12+EL22-4*EL33).*cos(a).^2+EL12-EL22+2*EL33);
EH22 = (cos(a).^4).*(EL11-2*EL12+EL22-4*EL33)+(cos(a).^2).*(-2*EL11+2*EL12+4*EL33)+EL11;
EH23 = -sin(a).*cos(a).*((EL11-2*EL12+EL22-4*EL33).*(cos(a).^2)-EL11+EL12+2*EL33);
EH33 = (cos(a).^4).*(-EL11+2*EL12-EL22+4*EL33)+(cos(a).^2).*(EL11-2*EL12+EL22-4*EL33)+EL33;
end

function [EH11,EH12,EH13,EH22,EH23,EH33,dEH11dMu1,dEH12dMu1,dEH13dMu1,dEH22dMu1,dEH23dMu1,dEH33dMu1,dEH11dMu2,dEH12dMu2,dEH13dMu2,dEH22dMu2,dEH23dMu2,dEH33dMu2] = getElasticityDerivativesRank2(a,muPhys,E,nu,EMin)
% Properties in local frame
denominator = (1- muPhys(:,2)+muPhys(:,1).*muPhys(:,2)*(1-nu^2));
EL11 = E./denominator.*(muPhys(:,1))+EMin/(1-nu^2);
EL12 = E./denominator.*(muPhys(:,2).*muPhys(:,1)*nu)+nu*EMin/(1-nu^2);
EL22 = E./denominator.*(muPhys(:,2).*(1-muPhys(:,2)+muPhys(:,2).*muPhys(:,1)))+EMin/(1-nu^2);
EL33 = EMin/(1-nu^2)*(1-nu)/2;
dEL11dMu1 = (E./(denominator.^2)).*(1-muPhys(:,2));
dEL12dMu1 = (E./(denominator.^2)).*muPhys(:,2)*nu.*(1-muPhys(:,2));
dEL22dMu1 = (E./(denominator.^2)).*(muPhys(:,2).^2*nu^2.*(1-muPhys(:,2)));
dEL33dMu1 = zeros(size(dEL11dMu1));
dEL11dMu2 = (E./(denominator.^2)).*(muPhys(:,1).*(1-muPhys(:,1)*(1-nu^2)));
dEL12dMu2 = (E./(denominator.^2)).*muPhys(:,1)*nu;
dEL22dMu2 = (E./(denominator.^2)).*(1-2*muPhys(:,2)+2*muPhys(:,1).*muPhys(:,2)+(muPhys(:,2)).^2-2*(muPhys(:,2)).^2.*muPhys(:,1)+(muPhys(:,1)).^2.*(muPhys(:,2)).^2+(muPhys(:,2)).^2*nu^2.*muPhys(:,1)-(muPhys(:,2)).^2.*(muPhys(:,1)).^2*nu^2);
dEL33dMu2 = zeros(size(dEL11dMu2));
% Rotate to global frame
[EH11,EH12,EH13,EH22,EH23,EH33] = rotateLocaToGlobal(a,EL11,EL12,EL22,EL33);
[dEH11dMu1,dEH12dMu1,dEH13dMu1,dEH22dMu1,dEH23dMu1,dEH33dMu1] = rotateLocaToGlobal(a,dEL11dMu1,dEL12dMu1,dEL22dMu1,dEL33dMu1);
[dEH11dMu2,dEH12dMu2,dEH13dMu2,dEH22dMu2,dEH23dMu2,dEH33dMu2] = rotateLocaToGlobal(a,dEL11dMu2,dEL12dMu2,dEL22dMu2,dEL33dMu2);
end

function []= plotFinalDesign(nelX,nelY,rho,muPhys,a)
rho = reshape(rho,nelY,nelX); 
muPhys1 = reshape(muPhys(:,1),nelY,nelX);
muPhys2 = reshape(muPhys(:,2),nelY,nelX);
muMax = 0.5*(max(max(abs(muPhys1(:)),abs(muPhys2(:))))*2);
a = -reshape(a,nelY,nelX); % Negative because y is pointing down in imagesc
colormap(gray); imagesc(-rho); axis equal; axis tight; axis off;pause(1e-6);
hold on
colormap(gray); quiver(cos(a).*muPhys1/muMax,sin(a).*muPhys1/muMax,0,'r.','LineWidth',4); axis equal; axis tight; axis off;
colormap(gray); quiver(-cos(a).*muPhys1/muMax,-sin(a).*muPhys1/muMax,0,'r.','LineWidth',4); axis equal; axis tight; axis off;drawnow;
colormap(gray); quiver(cos(a+pi/2).*muPhys2/muMax,sin(a+pi/2).*muPhys2/muMax,0,'b.','LineWidth',4); axis equal; axis tight; axis off;
colormap(gray); quiver(-cos(a+pi/2).*muPhys2/muMax,-sin(a+pi/2).*muPhys2/muMax,0,'b.','LineWidth',4); axis equal; axis tight; axis off;drawnow;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is intended for educational purposes and is supplementary to    %
% the review paper                                                         %
% Topology Optimization of Multi-scale Structures: A Review,               %
% Jun Wu, Ole Sigmund, Jeroen P. Groen                                     %
% Submitted to Structural and Multidisciplinary Optimization, 2021         %
%                                                                          %
% This code is based on                                                    %
% Efficient topology optimization in MATLAB using 88 lines of code,        %
% by E. Andreassen et al., Struct Multidisc Optim (2011).                  %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserve all rights but do not guaranty that the code is      %
% free from errors. Furthermore, we shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
