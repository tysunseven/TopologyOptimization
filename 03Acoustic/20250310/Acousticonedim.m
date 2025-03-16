% 参数设置
k = 2;           % 波数
lambda = 2*pi/k; % 波长
p_in = 1;        % 振幅
N = 450;         % 单元数目
L = 4.5;           % 区间长度
h = L/N;         % 单元长度
m = lambda/h;    % 一个波长中单元的个数
right_absorbing = true;

% 生成网格
x_nodes = linspace(0,L,N+1)';

% 初始化刚度矩阵和质量矩阵
K = zeros(N+1,N+1);
M = zeros(N+1,N+1);

% 组装刚度矩阵和质量矩阵
for e = 1:N
    nodes = [e,e+1];
    local_K = [1,-1;-1,1]/h;
    local_M = h*[2,1;1,2]/6;
    K(nodes,nodes) = K(nodes,nodes) + local_K;
    M(nodes,nodes) = M(nodes,nodes) + local_M;
end

% 构建系统矩阵和右端向量(左入射右反射)
A = -K+(k^2)*M;
A(1,1) = A(1,1)-1i*k;

if right_absorbing
    A(N+1,N+1) = A(N+1,N+1) - 1i*k; % 右端吸收条件
end

F = zeros(N+1,1);
F(1) = -2i*k*p_in;

% 求解线性方程组
P = A\F;

% 计算解析解
if right_absorbing
    p_analytical = p_in*exp(-1i*k*x_nodes);
else
    p_analytical = p_in*exp(-2i*k*L)*exp(1i*k*x_nodes) + p_in*exp(-1i*k*x_nodes);
end

p_analytical = real(p_analytical);

% 绘制结果
figure;
plot(x_nodes, real(P), 'b-', x_nodes, p_analytical, 'r--');
legend('数值解', '解析解', 'Location', 'Best');
xlabel('x');
ylabel('Re(p)');
title(['数值解与解析解比较 (L=', num2str(L), ')']);
grid on;

% 计算最大误差
error = max(abs(real(P) - p_analytical));
disp(['最大误差：', num2str(error)]);