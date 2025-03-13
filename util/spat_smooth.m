function out = spat_smooth(rate_map, sigma, type)
% 
if type == 1

    if sigma ~= 0
        [x, y] = meshgrid(-1:1, -1:1);
        filter = exp(-(x.^2 + y.^2) / (2*sigma^2)) / (2*pi*sigma^2);
        filter = filter / sum(filter(:));
        out = conv2(rate_map, filter);
    
        out = out(2:end-1, :);
        out = out(:, 2:end-1);
    else
        out = rate_map;
    end

else
    out = imgaussfilt(rate_map,sigma);
end