clearvars -global; close all; clc;
for j = 3:3
timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
target_mode = j; num_random_points = 0; num_examples = 30; nelx = 32;

% % 如果是随机生成初始
% output_folder = sprintf('x_init_data_mode%d_rand%d_nelx%d_%s', target_mode, num_random_points, nelx, timestamp); mkdir(output_folder);

% 如果是读取之前已有的初始
input_folder = 'x_init_data_mode1_rand0_nelx32_20250505_224641';
% output_folder = [input_folder '_retry_rand%d_%s',num_random_points, timestamp]; mkdir(output_folder);
output_folder = sprintf('x_retry_data_mode%d_rand%d_nelx%d_%s', target_mode, num_random_points, nelx, timestamp); mkdir(output_folder);

resultTable = table('Size',[num_examples 19],...
    'VariableTypes', [repmat({'double'},1,3),repmat({'double'},1,7), {'cell'}, repmat({'double'},1,7),{'double'}],...
    'VariableNames',{'loop','time','avgtime','bdmax','bdmin','bdmax_x','bdmax_y','bdmin_x','bdmin_y','bdgap',...
                    'EmptyColumn','gbmax','gbmin','gbmax_x','gbmax_y','gbmin_x','gbmin_y','gbgap','IsBoundaryExtrema'});

for i = 1:num_examples
    clearvars -except output_folder resultTable i timestamp target_mode num_random_points input_folder nelx

    % % 如果是随机生成初始
    % x_init = rand(nelx,nelx); x_init_filename = sprintf('rand%02d.xlsx', i); writematrix(x_init, fullfile(output_folder, x_init_filename));
    
    % 如果是读取之前已有的初始
    x_init_filename = sprintf('rand%02d.xlsx', i); x_path = fullfile(input_folder, x_init_filename); x_init = readmatrix(x_path); writematrix(x_init, fullfile(output_folder, x_init_filename));

    [x_final, figure, loop, time, bdmax, bdmin, bdmax_x, bdmax_y, bdmin_x, bdmin_y, bdgap, gbmax, gbmin, gbmax_x, gbmax_y, gbmin_x, gbmin_y, gbgap] = ...
        OutOfPlaneElasticBandInverseDesignMMaVideoMaxSqrt(x_init, target_mode, num_random_points, output_folder, i);

    x_final_filename = sprintf('rand%02d_final.xlsx', i); writematrix(x_final, fullfile(output_folder, x_final_filename));
    figure_filename = sprintf('rand%02d_figure.png', i); full_path = fullfile(output_folder, figure_filename); saveas(figure, full_path); 

    is_boundary_extrema = ~((abs(bdmax_x - gbmax_x) < 1e-6) && (abs(bdmax_y - gbmax_y) < 1e-6) && (abs(bdmin_x - gbmin_x) < 1e-6) && (abs(bdmin_y - gbmin_y) < 1e-6));

    resultTable(i,:) = {loop, time, time/loop, bdmax, bdmin, bdmax_x, bdmax_y, bdmin_x, bdmin_y, bdgap,...
                       {''}, gbmax, gbmin, gbmax_x, gbmax_y, gbmin_x, gbmin_y, gbgap, is_boundary_extrema};
end

% 如果是随机生成初始
result_filename = sprintf('FinalResult_retry_mode%d_rand%d_nelx%d_%s.xlsx', target_mode, num_random_points, nelx, timestamp);
writetable(resultTable, result_filename, 'WriteVariableNames', true, 'WriteMode', 'overwrite', 'Sheet', 'Results', 'Range', 'A1');

end