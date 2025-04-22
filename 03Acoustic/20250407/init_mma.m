function [m, n, xmin, xmax, xold1, xold2, low, upp, a0, a, c_MMA, d, move] = init_mma(nelx, nely, x)
% 优化参数初始化函数
% 输入参数:
%   nelx - X方向单元数
%   nely - Y方向单元数
%   x    - 初始设计变量矩阵
% 返回参数:
%   m, n, xmin, xmax, xold1, xold2, low, upp, a0, a, c_MMA, d, move

m     = 1;                % The number of general constraints. Volumn constrain
n     = nelx*nely;        % The number of design variables x_j.
xmin  = zeros( n, 1 );    % Column vector with the lower bounds for the variables x_j.
xmax  = ones( n, 1 );     % Column vector with the upper bounds for the variables x_j.
xold1 = x(:);          % xval, one iteration ago (provided that iter>1).
xold2 = x(:);          % xval, two iterations ago (provided that iter>2).
low   = [  ];             % Column vector with the lower asymptotes from the previous iteration (provided that iter>1).
upp   = [  ];             % Column vector with the upper asymptotes from the previous iteration (provided that iter>1).
a0    = 1;                % The constants a_0 in the term a_0*z.
a     = zeros(m,1);       % Column vector with the constants a_i in the terms a_i*z.
c_MMA = 1000*ones(m,1);   % Column vector with the constants c_i in the terms c_i*y_i.
d     = zeros(m,1);       % Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
move  = 0.2;              % Max allowed change in design variables
end