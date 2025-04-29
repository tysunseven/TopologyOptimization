clearvars -global
close all; clc;


target_mode = 1; num_random_points = 0;


timestamp = datestr(now, 'yyyymmdd_HHMMSS');

output_folder = sprintf('x_init_data_mode%d_rand%d_%s',...
                       target_mode, num_random_points, timestamp);
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

resultTable = table('Size',[30 19],...
    'VariableTypes', [repmat({'double'},1,3),repmat({'double'},1,7), {'cell'}, repmat({'double'},1,7),{'double'}],...
    'VariableNames',{'loop','time','avgtime','bdmax','bdmin','bdmax_x','bdmax_y','bdmin_x','bdmin_y','bdgap',...
                    'EmptyColumn',...  % 空列占位
                    'gbmax','gbmin','gbmax_x','gbmax_y','gbmin_x','gbmin_y','gbgap','IsBoundaryExtrema'});

for i = 1:5
    % 生成新随机初始设计变量
    clearvars -except output_folder resultTable i timestamp target_mode num_random_points input_folder

    x_init = rand(30,30);
    x_filename = sprintf('rand%02d.xlsx', i);
    writematrix(x_init, fullfile(output_folder, x_filename));

    % 运行优化算法
    [x_final, figure, loop, time, bdmax, bdmin, bdmax_x, bdmax_y, bdmin_x, bdmin_y, bdgap, gbmax, gbmin, gbmax_x, gbmax_y, gbmin_x, gbmin_y, gbgap] = ...
        OutOfPlaneElasticBandInverseDesignMMaVideoMaxSqrt(x_init, target_mode, num_random_points);

    x_filename = sprintf('rand%02d_final.xlsx', i);
    writematrix(x_final, fullfile(output_folder, x_filename));

    figure_filename = sprintf('rand%02d_figure.png', i);
    full_path = fullfile(output_folder, figure_filename);
    saveas(figure, full_path); 

    tol = 1e-6; 

    % 检查最大值位置是否相同
    max_pos_match = (abs(bdmax_x - gbmax_x) < tol) && ...
                   (abs(bdmax_y - gbmax_y) < tol);

    % 检查最小值位置是否相同
    min_pos_match = (abs(bdmin_x - gbmin_x) < tol) && ...
                   (abs(bdmin_y - gbmin_y) < tol);

    % 综合判断标识 (0:完全匹配, 1:存在不匹配)
    is_boundary_extrema = ~(max_pos_match && min_pos_match);



    resultTable(i,:) = {loop, time, time/loop, bdmax, bdmin, bdmax_x, bdmax_y, bdmin_x, bdmin_y, bdgap,...
                       {''},...  % 空列填充空字符串
                       gbmax, gbmin, gbmax_x, gbmax_y, gbmin_x, gbmin_y, gbgap, is_boundary_extrema};
end


result_filename = sprintf('FinalResult_mode%d_rand%d_%s_retry2_rand30.xlsx',...
                          target_mode, num_random_points, timestamp);
writetable(resultTable, result_filename,...
          'WriteVariableNames', true,...
          'WriteMode', 'overwrite',...
          'Sheet', 'Results',...
          'Range', 'A1');