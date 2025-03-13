function [out, foo] = state_mean(state, inp)
% computes nanmean+std+num of inp for each unique state in 'state'
% input: STATE: vector of state label [1 1 2 3 2 1]
%          INP: vector of corresponding data [0.1 0.2 0.3 0.1 0.2 0.1] 
% output: OUT: matrix of [state, mean, num, std] 
%         FOO: nan rows in inp which were removed
% automatically determines if 'inp' has a time component (assumes 2nd dim is time)

%      keyboard

    out = [];
    indx = 1;

    if size(state, 1) == 1
        ... % don't need to change it
    elseif size(state, 2) == 1
        state = state';
    end

    if size(inp, 1) == 1 && size(inp, 2) > 1
        inp = inp';
    end

    if size(inp, 3) == 1

        foo = find(any(sum(isnan(inp), 2) > 1));
        rel = setdiff(1:size(inp, 1), foo);

        state = state(rel);
        if size(inp, 2) > 1
           inp   = inp(rel, :);
        else
           inp   = inp(rel);
        end
    elseif size(inp, 3) > 1
        foo = find(any(sum(isnan(inp), 2) > 1, 3));
        rel = setdiff(1:size(inp, 1), foo);

        state = state(rel);
        inp   = inp(rel, :, :);
    end

    if size(inp, 2) > 1
         for e = unique(state)
            rel_st  = e == state;
            tmp     = nanmean(inp(find(rel_st), :, :), 1);
            tmp_std = nanstd(inp(find(rel_st), :, :), [], 1) ./ (sqrt(sum(rel_st)));
            if size(inp, 3) > 1 % we're doing this on a 3d mat
                out(indx, 1, :, :) = repmat(e, [1, size(tmp, 2), size(tmp, 3)]);
                out(indx, 2, :, :) = tmp;
                out(indx, 3, :, :) = repmat(sum(rel_st), [1, size(tmp, 2), size(tmp, 3)]);
                out(indx, 4, :, :) = tmp_std;
            else % we're doing this on a [trials x time matrix]
                out(indx, 1, :) = repmat(e, [1, size(tmp, 2), size(tmp, 3)]);
                out(indx, 2, :) = tmp;
                out(indx, 3, :) = repmat(sum(rel_st), [1, size(tmp, 2), size(tmp, 3)]);
                out(indx, 4, :) = tmp_std;
            end
            indx = indx + 1;
         end
    else
        for e = unique(state)
            rel_st  = e == state;
            tmp     = nanmean(inp(find(rel_st)));
            tmp_std = nanstd(inp(find(rel_st))) / (sqrt(sum(rel_st)));
            out(indx, :) = [e, tmp, sum(rel_st), tmp_std];
            indx = indx + 1;
        end
    end

    if size(inp, 2) > 1
        % keyboard
        rm_rw = find(~isnan(out(:, 1, 1, 1)));
        out = out(rm_rw, :, :, :);
    else
        rm_rw = find(~isnan(out(:, 1)));
        out = out(rm_rw, :);
    end

end
