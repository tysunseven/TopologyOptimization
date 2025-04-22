function [row, col, fixT] = init_trans(nelx, nely)

row = [1,((nely+1)+1):(nely+1):((nelx-1)*(nely+1)+1),nelx*(nely+1)+1,(nelx*(nely+1)+2):((nelx+1)*(nely+1))];
col = [nely+1,(2*(nely+1):(nely+1):(nelx*(nely+1))),nely+1,2:nely,nely+1];
fixdofs = setdiff(1:(nely+1)*(nelx+1),row);
fixT = sparse(fixdofs,fixdofs,ones(size(fixdofs)),(nely+1)*(nelx+1),(nely+1)*(nelx+1));

end

