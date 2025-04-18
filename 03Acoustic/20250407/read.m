x = imread("g6.jpg");
x = x(:,:,1);
x = x/255;
x = 1-x;
disp(x);
% valid_rows = any(x, 2); % 找出包含非零元素的行
% valid_cols = any(x, 1); % 找出包含非零元素的列
% x = x(valid_rows, valid_cols);
% [m, n] = size(x);
% max_dim = max(m, n);               % 取行列最大值
% pad_x = max_dim - m;               % 需要补的行数
% pad_y = max_dim - n;               % 需要补的列数
% 
% % 创建方阵并填充（补在右侧和下侧）
% x_square = zeros(max_dim, max_dim, 'like', x);
% x_square(1:m, 1:n) = x;
% x=x_square;
disp(size(x));

figure('Position', [100 100 400 400]);
colormap(gray);
imagesc(1-x);
clim([0 1]);
axis equal tight;
axis off;
title('Material Distribution');
drawnow;
