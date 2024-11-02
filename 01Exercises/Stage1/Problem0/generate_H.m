function H = generate_H(rmin)
    % Parameters
    nelx = 3;
    nely = 3;
    
    % Precompute the total number of elements
    % iH, jH, and sH are initialized according to the given formula
    num_elements = nelx * nely * (2*(ceil(rmin)-1)+1)^2;
    iH = ones(num_elements, 1);
    jH = ones(num_elements, 1);
    sH = zeros(num_elements, 1);
    
    % Counter for filling iH, jH, and sH
    k = 0;
    
    % Nested loops to calculate the values in H
    for i1 = 1:nelx
        for j1 = 1:nely
            e1 = (i1-1) * nely + j1;
            for i2 = max(i1 - (ceil(rmin)-1), 1):min(i1 + (ceil(rmin)-1), nelx)
                for j2 = max(j1 - (ceil(rmin)-1), 1):min(j1 + (ceil(rmin)-1), nely)
                    e2 = (i2-1) * nely + j2;
                    k = k + 1;
                    iH(k) = e1;
                    jH(k) = e2;
                    sH(k) = max(0, rmin - sqrt((i1-i2)^2 + (j1-j2)^2));
                end
            end
        end
    end
    
    % Construct the sparse matrix H
    H = sparse(iH, jH, sH);
    Hs = sum(H,2);
    HsH = spdiags(1.0./Hs,0,nelx*nely,nelx*nely)*H;
    A = (HsH'*repmat(1/(nelx*nely),nelx*nely,1))';
    A = reshape(A, 3, 3);
    dv = ones(nely,nelx);
    % H的每一行是对应单元的权重信息
    % Hs是行的和
    % dv是按列展平
    % Output the matrix H
    disp('Matrix H:');
    full(H) % Display the full matrix for better visualization
    disp('Hs:');
    disp(Hs);
    disp('dv(:)./Hs:');
    disp(dv(:)./Hs);
    disp('H*(dv(:)./Hs)');
    disp(H*(dv(:)./Hs));
    disp('Matrix HsH:');
    full(HsH) % Display the full matrix for better visualization
    disp('A:');
    full(A) % Display the full matrix for better visualization
    dv(:)=H*(dv(:)./Hs);
    dv=reshape(dv,3,3);
    dv = dv / 9;
    disp('dv:');
    disp(dv);
end