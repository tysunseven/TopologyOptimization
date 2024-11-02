% 定义区间和步长
x = 0:0.01:1;

% 不同的 beta 值
beta_values = [1, 2, 4, 8, 16, 32];

figure;
hold on;

% 遍历不同的 beta 值并绘制函数图像
for i = 1:length(beta_values)
    beta = beta_values(i);
    y = exp(-beta * (1-x)) - (1-x) .* exp(-beta);
    plot(x, y, 'LineWidth', 1.5);
end

legend('\beta = 1', '\beta = 2', '\beta = 4', '\beta = 8', '\beta = 16', '\beta = 32', 'Location', 'NorthWest');
xlabel('x');
ylabel('y');
title('Plot of y = e^{-\beta (1-x)} - (1-x) e^{-\beta}');
grid on;
hold off;

% 保存为 PNG 图像文件
saveas(gcf, 'Heaviside2.png');