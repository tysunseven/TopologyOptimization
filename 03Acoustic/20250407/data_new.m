clearvars -global; close all; clc;
timestamp = datetime('now', 'Format', 'yyyyMMdd_HHmmss');

target_mode = 1; num_random_points = 0; number = 30;

output_folder = sprintf('x_init_data_mode%d_rand%d_%s',...
                       target_mode, num_random_points, timestamp);
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

resultTable = table('Size',[number 19],...
    'VariableTypes', [repmat({'double'},1,3),repmat({'double'},1,7), {'cell'}, repmat({'double'},1,7),{'double'}],...
    'VariableNames',{'loop','time','avgtime','bdmax','bdmin','bdmax_x','bdmax_y','bdmin_x','bdmin_y','bdgap',...
                    'EmptyColumn','gbmax','gbmin','gbmax_x','gbmax_y','gbmin_x','gbmin_y','gbgap','IsBoundaryExtrema'});

for i = 1:number

    clearvars -except output_folder resultTable i timestamp target_mode num_random_points number

    x_init = rand(30,30); x_filename = sprintf('rand%02d.xlsx', i);
    writematrix(x_init, fullfile(output_folder, x_filename));

    [x_final, figur, loop, time, bdmax, bdmin, bdmax_x, bdmax_y, bdmin_x, bdmin_y, bdgap, gbmax, gbmin, gbmax_x, gbmax_y, gbmin_x, gbmin_y, gbgap] = ...
        OutOfPlaneElasticBandInverseDesignMMaVideoMaxSqrt(x_init, target_mode, num_random_points, output_folder, i);

    x_filename = sprintf('rand%02d_final.xlsx', i);
    writematrix(x_final, fullfile(output_folder, x_filename));

    figure_filename = sprintf('rand%02d_figure.png', i);
    full_path = fullfile(output_folder, figure_filename);
    saveas(figur, full_path); 

    x_repmat = repmat(x_final, 3, 3);
    figure;  % <- 关键修改：隐藏图形窗口
    imagesc(1-x_repmat);                     % 显示矩阵颜色图
    colormap(gray); axis off;
    figure_filename = sprintf('rand%02d_repmat.png', i);
    full_path = fullfile(output_folder, figure_filename);
    saveas(gcf, full_path); 

    tol = 1e-6; 

    max_pos_match = (abs(bdmax_x - gbmax_x) < tol) && ...
                   (abs(bdmax_y - gbmax_y) < tol);
    min_pos_match = (abs(bdmin_x - gbmin_x) < tol) && ...
                   (abs(bdmin_y - gbmin_y) < tol);
    is_boundary_extrema = ~(max_pos_match && min_pos_match);

    resultTable(i,:) = {loop, time, time/loop, bdmax, bdmin, bdmax_x, bdmax_y, bdmin_x, bdmin_y, bdgap,...
                       {''}, gbmax, gbmin, gbmax_x, gbmax_y, gbmin_x, gbmin_y, gbgap, is_boundary_extrema};
end


result_filename = sprintf('FinalResult_mode%d_rand%d_%s.xlsx',...
                          target_mode, num_random_points, timestamp);
writetable(resultTable, result_filename,...
          'WriteVariableNames', true,...
          'WriteMode', 'overwrite',...
          'Sheet', 'Results',...
          'Range', 'A1');