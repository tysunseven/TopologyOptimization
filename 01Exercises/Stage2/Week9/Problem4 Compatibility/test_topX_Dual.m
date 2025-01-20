% Set parameters
volfrac = 0.5; % Set the volume fraction
nely = 100; % Set the number of elements in the y-direction
nelx = 100; % Set the number of elements in the x-direction

% Create a figure
figure;

for IC = 1:3
    % Initialize design variable matrix
    x = repmat(volfrac, nely, nelx);
    
    % Generate initial condition pattern
    for i = 1:nelx
        for j = 1:nely
            if IC == 1
                if sqrt((i - nelx/2 - 0.5)^2 + (j - nely/2 - 0.5)^2) < min(nelx, nely) / 6
                    x(j, i) = volfrac / 2;
                end
            elseif IC == 2
                if sqrt((i - nelx/4 - 0.5)^2 + (j - nely/4 - 0.5)^2) < min(nelx, nely) / 6
                    x(j, i) = volfrac / 2;
                end
            elseif IC == 3
                if sqrt((i - nelx/4 - 0.5)^3 + (j - nely/8 - 0.5)^2) < min(nelx / 2, nely) / 6
                    x(j, i) = volfrac / 2;
                end
            end
        end
    end

    % Plot each pattern in a subplot
    subplot(1, 3, IC);
    imagesc(1 - x, [0, 1]);
    colormap(gray);
    title(['IC = ', num2str(IC)]);
    axis equal; % 保持纵横比相等
    xlim([0, nelx]); % 限制x轴显示范围
    ylim([0, nely]); % 限制y轴显示范围
    axis off; % 隐藏坐标轴
end