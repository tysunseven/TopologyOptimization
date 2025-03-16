% 固定参数
L = 5;                  % 固定区间长度
right_absorbing = true; % 固定吸收边界条件
p_in = 1;               % 固定振幅

% 参数扫描范围
k_values = [0.5, 1, 2, 4, 8];      % 变化的波数
N_values = [50, 100, 200, 400, 800]; % 变化的网格数

% 预分配存储矩阵
results = zeros(length(k_values)*length(N_values), 4);
row = 1;

% 参数扫描主循环
for ki = 1:length(k_values)
    k = k_values(ki);
    
    for Ni = 1:length(N_values)
        N = N_values(Ni);
        
        % 计算相关参数
        h = L/N;
        lambda = 2*pi/k;
        m = lambda/h;
        
        % 运行仿真
        [error, ~] = helmholtz_solver(k, p_in, N, L, right_absorbing);
        
        % 存储结果
        results(row,:) = [k, N, m, error];
        row = row + 1;
    end
end

% 创建结果表格并排序
result_table = array2table(results, ...
    'VariableNames', {'k','N','m','error'});
result_table = sortrows(result_table, 'N');  % 按N从小到大排序

% 显示格式化表格
disp('=================== 参数扫描结果 (按N排序) ===================');
disp('    k       N            m        最大误差');
disp('-----------------------------------------------------------');
for i = 1:size(result_table,1)
    fprintf('%6.2f  %6d  %12.2f  %12.4e\n',...
            result_table.k(i),...
            result_table.N(i),...
            result_table.m(i),...
            result_table.error(i));
end

%% 求解器函数保持不变
function [error, P] = helmholtz_solver(k, p_in, N, L, right_absorbing)
    h = L/N;
    x_nodes = linspace(0,L,N+1)';

    K = sparse(N+1,N+1);
    M = sparse(N+1,N+1);
    
    for e = 1:N
        nodes = [e,e+1];
        K(nodes,nodes) = K(nodes,nodes) + [1,-1;-1,1]/h;
        M(nodes,nodes) = M(nodes,nodes) + h*[2,1;1,2]/6;
    end

    A = -K + (k^2)*M;
    A(1,1) = A(1,1) - 1i*k;
    
    if right_absorbing
        A(N+1,N+1) = A(N+1,N+1) - 1i*k;
    end

    F = sparse(N+1,1);
    F(1) = -2i*k*p_in;
    P = A\F;

    p_analytical = p_in*exp(-1i*k*x_nodes);
    p_analytical = real(p_analytical);

    error = max(abs(real(P) - p_analytical));
end