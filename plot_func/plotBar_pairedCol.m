function out = plotBar_pairedCol(inp_val, col_a, col_b, version)

for tp = 1:2

    if tp == 1
        rel_vals = inp_val(:, 1:2);
        rel_x_base = 0;
    else
        rel_vals = inp_val(:, 3:4);
        rel_x_base = 2;
    end

    % plot dots
    for n_p = 1:size(rel_vals, 1)
        vals = [];
        X = find(~isnan(rel_vals(n_p, :)));
        
        for e = 1:length(X)
            vals(e) = rand * 0.15 + 0.025;
        end

        X = X + vals + rel_x_base;
        Y = rel_vals(n_p, :);
        Y = Y(~isnan(Y));
        alph_lev = 0.7;
        line(X, Y, 'Color', [0, 0, 0, 0.3], 'LineWidth', 1); hold on

        for e = 1:size(X, 2)

            if tp == 1
                rel_curr_col_a = col_a;
                if version == 1
                    rel_curr_col_b = col_a;
                else
                    rel_curr_col_b = col_a - 0.25;
                end
            else
                rel_curr_col_a = col_b;
                if version == 1
                    rel_curr_col_b = col_b;
                else
                    rel_curr_col_b = col_b - 0.25;
                end
            end


            scatter(X(e), Y(e), 'filled', ... 
                      'MarkerEdgeColor', [0.2 0.2 0.2],     ...
                      'MarkerFaceColor', rel_curr_col_a, ...
                      'MarkerFaceAlpha',  alph_lev, ...
                      'MarkerEdgeAlpha', alph_lev, ... 
                      'LineWidth', 1); hold on
        end
    end

    nStates = size(rel_vals, 2);
    b = bar([1, 2] + rel_x_base, nanmean(rel_vals, 1), 'edgeColor', 'k', 'LineWidth', 3); hold on
    b.FaceColor = 'flat';
    
    % we select which one we want to color
    b.CData(1, :) = repmat(rel_curr_col_a, [1, 1]);
    b.CData(2, :) = repmat(rel_curr_col_b, [1, 1]);
    err = nanstd(rel_vals, [], 1)./sqrt(sum(~isnan(rel_vals), 1));
    errorbar([1, 2] + rel_x_base, nanmean(rel_vals, 1), err,  '.', 'CapSize', 0, 'linewidth', 3, 'color', [1/255, 1/255, 1/255])
    b.FaceAlpha = 0.7;
end