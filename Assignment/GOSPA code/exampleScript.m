% AUTHOR: Abu Sajana Rahmathullah
% DATE OF CREATION: 7 August, 2017

% set parameters
p = 2; c = 8; alpha = 2; ndim = 2;

% Generate random sets
% columns of x_mat and y_mat correspond to different points in the sets
nx      = randi([0, 10]);
ny      = randi([0, 10]);
x_mat   = 20 * (rand(ndim, nx) - 0.5);
y_mat   = 20 * (rand(ndim, ny) - 0.5);

% compute GOSPA
[d_gospa, x_to_y_assignment, decomp_cost] = GOSPA(x_mat, y_mat, p, c, alpha);

% plot the input and the output if the vectors are two dimensional
if ndim == 2
    % plot the input vectors
    x0 = 0; y0 = 0; width = 10; height = 4;
    f1 = figure('Units','inches','Position',[x0 y0 width height], ...
        'PaperPositionMode','auto');
    hold on;
    bluecol = [0.3 0.3 1];
    redcol  = [1 0.2 0.2];
    text(x_mat(1,:), x_mat(2,:), num2str((1:nx)'), ...
        'HorizontalAlignment', 'right');
    text(y_mat(1,:), y_mat(2,:), num2str((1:ny)'), ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'bottom');
    hx = plot(x_mat(1,:), x_mat(2,:), 'x', 'color', bluecol, 'MarkerSize', 10);
    hy = plot(y_mat(1,:), y_mat(2,:), 'o', 'color', redcol, 'MarkerSize', 10);
    
    % plot the output
    x_ind = 1:nx;
    x_ind = x_ind(x_to_y_assignment~=0);
    y_ind = x_to_y_assignment(x_to_y_assignment~=0);
    for ind = 1:sum(x_to_y_assignment~=0)
        db = min(sqrt(sum((x_mat(:,x_ind(ind)) - y_mat(:,y_ind(ind))).^2)),c);        
        if db < c
            plot([x_mat(1, x_ind(ind)), y_mat(1, y_ind(ind))], ...
                 [x_mat(2, x_ind(ind)), y_mat(2, y_ind(ind))], 'k-');
            
            text((x_mat(1, x_ind(ind))+y_mat(1, y_ind(ind)))/2, ...
                 (x_mat(2, x_ind(ind))+y_mat(2, y_ind(ind)))/2, ...
                 num2str(db, 3), 'HorizontalAlignment', 'right', ...
                 'VerticalAlignment', 'bottom');
        end
    end
    hlegend = legend([hx, hy], {'X', 'Y'});
    set(hlegend, 'Location', 'Best');
    set(gca, 'FontSize', 18, 'FontName', 'Times');
    title(['c=' num2str(c), ', p=' num2str(p), ', \alpha=' ...
        num2str(alpha) ', GOSPA=' num2str(d_gospa, 3), ]);
end