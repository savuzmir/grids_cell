function out = imagesc_text(inp, inp_txt)

%   keyboard
    imagesc(inp);

    imAlpha=ones(size(inp));
    imAlpha(isnan(inp))=0;
    imagesc(inp,'AlphaData',imAlpha);
    set(gca,'color',1*[1 1 1]);

    N = size(inp, 1);
    x = repmat(1:N,N,1); y = x';
    t = num2cell(round(inp_txt*100, 0)); % extact values into cells
    t = cellfun(@num2str, t, 'UniformOutput', false); % convert to string

    col_corr = [255, 255, 255] ./ 255;% [166, 184, 244, 245]./255;
    col_corr = repmat(col_corr, [length(x(:)), 1]);
    col_corr(find(imAlpha(:) == 0), :) = repmat([1, 1, 1, 1], [length(find(imAlpha(:) == 0)), 1]);

    tmp_x = x(:);
    tmp_y = y(:);

    for xo = 1:length(x(:)) 
        text(x(xo), y(xo), t(xo), 'HorizontalAlignment', 'Center', 'color', col_corr(xo, :), 'FontSize', 14);
    end

end