function [K, M] = loop_fem(mu1,mu2,rho1,rho2,x,iIndex,jIndex,Ke,Me,nelx,nely)

mu = mu1*(1-x)+mu2*x;
rho = rho1*(1-x)+rho2*x;
sK = reshape(Ke(:)*(mu(:))',16*nelx*nely,1);
sM = reshape(Me(:)*(rho(:))',16*nelx*nely,1);
K = sparse(iIndex,jIndex,sK); K = (K+K')/2;
M = sparse(iIndex,jIndex,sM); M = (M+M')/2;

end