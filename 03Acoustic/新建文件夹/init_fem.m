function [edofMata, Ke, Me, iIndex, jIndex] = init_fem(nelx, nely,h)

nodenrsa = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVeca = reshape(nodenrsa(1:end-1,1:end-1)+nely+1,nelx*nely,1);
edofMata = repmat(edofVeca,1,4)+repmat([0 -nely-1 -nely 1],nelx*nely,1);

Ke = [4 -1 -2 -1; -1 4 -1 -2; -2 -1 4 -1; -1 -2 -1 4]/6;
Me = h^2*[4 2 1 2; 2 4 2 1; 1 2 4 2; 2 1 2 4]/36;

iIndex = reshape(kron(edofMata,ones(4,1))',16*nelx*nely,1);
jIndex = reshape(kron(edofMata,ones(1,4))',16*nelx*nely,1);

end