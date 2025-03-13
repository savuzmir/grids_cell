function plot_sigline(curr, pvals, color, position)

figure(curr)
% weight according to sig
b1_05  = nan(length(pvals), 1);
b1_01  = nan(length(pvals), 1);
b1_001 = nan(length(pvals), 1);

b1_05(pvals <= 0.05) = position;
b1_01(pvals <= 0.01) = position;
b1_001(pvals<=0.001) = position;

%b1_05(pvals>0.01 & pvals <= 0.05) = position;
%b1_01(pvals>0.001 & pvals <= 0.01) = position;
%b1_001(pvals<0.001) = position;
    
%keyboard

plot(1:length(pvals), b1_05, 'color', color, 'linewidth', 1.5); hold on
plot(1:length(pvals), b1_01, 'color', color, 'linewidth', 2.5); hold on
plot(1:length(pvals), b1_001, 'color', color, 'linewidth', 3.5)

% scatter(find(pvals <= 0.05), b1_05(pvals <= 0.05), 5.5, 'square', 'MarkerEdgeColor', color, 'MarkerFaceColor', color);
% scatter(find(pvals <= 0.01), b1_01(pvals <= 0.01), 10.5, 'square', 'MarkerEdgeColor', color, 'MarkerFaceColor', color);
% scatter(find(pvals <= 0.001), b1_001(pvals <= 0.001), 20.5, 'square', 'MarkerEdgeColor', color, 'MarkerFaceColor', color);


end