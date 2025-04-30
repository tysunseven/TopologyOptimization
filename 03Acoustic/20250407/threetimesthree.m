clearvars -global
close all; clc;



% 原始数据路径
input_folder = 'round2/x_init_data_mode3_rand0_20250424_120046_retry2_rand30';

for i = 1:30 

    x_filename = sprintf('rand%02d_final.xlsx', i);
    x_path = fullfile(input_folder, x_filename);
    x_init = readmatrix(x_path);
    x = repmat(x_init, 3, 3);
    
    % 绘制矩阵图像
    fig = figure('Visible', 'off');  % <- 关键修改：隐藏图形窗口
    imagesc(1-x);                     % 显示矩阵颜色图
    colormap(gray);
    
    % 生成保存路径（与原Excel文件同名，扩展名改为png）
    [~, name, ~] = fileparts(x_filename);  % 提取文件名（不含扩展名）
    output_filename = [name, '.png'];       % 新文件名
    output_path = fullfile(input_folder, output_filename);
    
    % 保存图像
    saveas(gcf, output_path, 'png');        % 保存为PNG
    
    close;
end

