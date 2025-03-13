function out = plot_dots(input, curr, figInfo, col_A, paired_samples)

if paired_samples
    nStates = size(input(:, 1));
else
    ...
end

all_X = [];
colo  = [0.5, 0.5, 0.5];

for i = 1:nStates
    if paired_samples

        vals = [];
        X = find(~isnan(input(i, :)));

        X(:, 2:end) = X(:, 2:end) / 1.75;
        
        for e = 1:length(X)
            vals(e) = rand * 0.03 + 0.025;
        end

        X = X + vals;
        Y = input(i, :);
        Y = Y(~isnan(Y));
    else
        X = zeros(sum(~isnan(input(:, i))), 1);
        Y = input(:, i);
        Y = Y(~isnan(Y));
    end

    alph_lev = 0.3;
    
    if isempty(Y)
        continue
    end

    if paired_samples
        line(X, Y, 'Color', [colo, alph_lev], 'LineWidth', 0.8);
    end

    for e = 1:size(X, 2)
        scatter(X(e), Y(e), 'filled', ... 
                  'MarkerEdgeColor', [0.2, 0.2, 0.2],     ...
                  'MarkerFaceColor', col_A{e}, ...
                  'MarkerFaceAlpha',  alph_lev, ...
                  'MarkerEdgeAlpha', alph_lev, ... 
                  'LineWidth', 1.5); hold on
    end

    all_X(i, :) = X;

end

% plot mean + sem
x_pos_available = 1:size(input, 2);
x_pos_available(2:end) = x_pos_available(2:end) / 1.75;

for ip = 1:size(input, 2)

    x_pos   = x_pos_available; x_pos = x_pos(ip);
    st_mean = nanmean(input(:, ip), 1);
    st_sem  = sem(input(:, ip), 1);
    x_lan = [x_pos; x_pos];
    y_lan = [st_mean-st_sem, st_mean+st_sem];
    plot(x_lan, y_lan, 'k', 'linewidth', 2.5); hold on
    scatter(x_pos, st_mean, 150, 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col_A{ip}, 'LineWidth', 2.5); hold on

end

% plot total mean
tot_mean = nanmean(nanmean(input, 1), 2);

midpoints = [];
for g = 1:size(all_X, 2)
    plot([prctile(all_X(:, g), 97.5, 1), x_pos_available(:, g)], [tot_mean, tot_mean], 'linewidth', 1.5, 'linestyle', '--', 'color', 'k');
    midpoints(g) = nanmean([prctile(all_X(:, g), 97.5, 1), x_pos_available(:, g)]);
end

figInfo{6} = midpoints;
figElements(curr, figInfo{:})
box off

end
