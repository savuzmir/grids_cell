function out = compute_rate_map(fr, mag_field, prob_field, overlap, overlap_step, smoothing_kernel, smothing_kernel_type, varargin)
% function out = compute_rate_map(fr, mag_field, prob_field, overlap, overlap_step, smoothing_kernel, smothing_kernel_type, varargin)
% generates  the rate map and convolves the matrix with a 3x3 gaussian filter.
% nans + infs are zeroed out; throws a warning if this happens so its known. very high/low tvalues (e.g. -300/300 happening due to low trials are set to -3/3 for visual purposes
% returns out 
% [dim1 x dim2 rate map; dim1 increases from smallest to largest across ROWS
%                        dim2 increases from smallest to largest across COLUMNS]

    rel_mag   = mag_field';
    rel_prob  = prob_field';

     unq_lv_mag  = unique(rel_mag)';
     unq_lv_prob = unique(rel_prob)';

    if ~isempty(varargin) % this argument is used for split data to ensure we always get equal coverage across sessions
        unq_lv_mag  = unique(varargin{1});
        unq_lv_prob = unique(varargin{2});
    end

    if length(unq_lv_mag) > 10 % we know it has to be value 

        if isempty(varargin)
            if size(unq_lv_mag, 1) > 1 & size(unq_lv_mag, 2) == 1
                unq_lv_mag  = unique([unq_lv_mag; unq_lv_prob]);
            else
                unq_lv_mag  = unique([unq_lv_mag, unq_lv_prob]);
            end
    
            unq_lv_prob = unq_lv_mag;
        else
            unq_lv_mag  = unique(varargin{1});
            unq_lv_prob = unique(varargin{2});
        end
    end

    curr_rate_map        = nan(length(unq_lv_mag) - overlap, length(unq_lv_mag) - overlap);
    curr_rate_stat_map   = curr_rate_map;
    curr_rate_map_sem    = curr_rate_map;
    curr_rate_map_counts = curr_rate_map; 

    indx_x = 1;

    for mg = 1:overlap_step:(length(unq_lv_mag) - overlap)
        rel_mg = unq_lv_mag(mg):unq_lv_mag(mg+overlap);
        indx_y  = 1;
        for pr = 1:overlap_step:(length(unq_lv_prob) - overlap)
            rel_pr = unq_lv_prob(pr):unq_lv_prob(pr+overlap);

            rel_mg_rows = find_rows(rel_mag', rel_mg, 'yes');
            rel_pr_rows = find_rows(rel_prob', rel_pr, 'yes');

            if iscell(rel_mg_rows)
                joint_rows  = intersect(rel_mg_rows{1}, rel_pr_rows{1});
            else
                joint_rows = intersect(rel_mg_rows, rel_pr_rows);
            end

            foo         = fr(joint_rows);
            foo         = foo(~isnan(foo));
            num_samples = nanlength(foo);

            curr_rate_map(indx_x, indx_y)        = nanmean(foo);
            [~,~,~,tstats]                       = ttest(foo, repmat(nanmean(fr), [length(foo), 1])); % against pop average resp

            try
                curr_rate_stat_map(indx_x, indx_y)   = tstats.tstat;
            catch
                curr_rate_stat_map(indx_x, indx_y)   = nan;
            end
            curr_rate_map_sem(indx_x, indx_y)    = nanstd(foo)./sqrt(num_samples);
            curr_rate_map_counts(indx_x, indx_y) = num_samples;
            indx_y = indx_y + 1;
        end
        indx_x = indx_x + 1;
    end

%     keyboard

    out = struct;

    if any(find(isnan(curr_rate_map))) 
        warning('added nan-zeros to rate map')

        % remove nan 
        if find(isnan(curr_rate_map))
            out.added_nan                               = 1;
            out.added_nan_pos                           = find(isnan(curr_rate_map));
            curr_rate_map(find(isnan(curr_rate_map)))   = 0;
        end
    end
    
    if any(find(isinf(curr_rate_map))) 
        warning('added inf-zeros to rate map')

        % remove inf 
        if find(isinf(curr_rate_map))
            out.added_inf                               = 1;
            out.added_inf_pos                           = find(isinf(curr_rate_map));
            curr_rate_map(find(isinf(curr_rate_map)))   = 0;
        end
    end

    curr_rate_map      = spat_smooth(curr_rate_map, smoothing_kernel, smothing_kernel_type);
    curr_rate_map_sem  = spat_smooth(curr_rate_map_sem, smoothing_kernel, smothing_kernel_type);

    out.rate_map_avg       = curr_rate_map;
    out.rate_map_sem       = curr_rate_map_sem;
    out.rate_stat_map      = curr_rate_stat_map;    
    out.rate_map_num       = curr_rate_map_counts;
end