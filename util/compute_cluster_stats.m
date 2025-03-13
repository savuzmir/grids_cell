function out = compute_cluster_stats(trueMatrix, permutationMatrix, thresh, dim)
% compute cluster-based permutation info
% trueMatrix (matrix of true values), permutationMatrix (matrix of permuted values)
% thresh - significance threshold to use
% dim - dimension along which we want to compute it 

% note: clustMassLB, trueMassLB returned are inverted with a -1, this is only important relevant 
% for permutation testing where pos and neg values exist in the data (e.g. cells)

    % out is all the relevant outputs
    threshUB = prctile(permutationMatrix, thresh, dim);
    threshLB = prctile(permutationMatrix, 100-thresh, dim);

    % find clusters of ones that exceed the permutation threshold 
    aboveBounds = permutationMatrix > threshUB;
    belowBounds = permutationMatrix < threshLB;
    
    threshLenUB   = [];
    threshMassUB  = [];
    fullIndex     = 1;
    
    % find living islands of above
    for i = 1:size(permutationMatrix, dim)

        if dim == 2
            rolledOutPermutationMatrix = squeeze(permutationMatrix(:, i)); rolledOutPermutationMatrix = rolledOutPermutationMatrix(:);
            [pockets, numPockets]      = bwlabel(aboveBounds(:,i));
        elseif dim == 3
            rolledOutPermutationMatrix = squeeze(permutationMatrix(:, :, i)); rolledOutPermutationMatrix = rolledOutPermutationMatrix(:);
            [pockets, numPockets]      = bwlabel(aboveBounds(:,:,i));
        end

        for pock = 1:numPockets
            tmpLen(fullIndex)  = length(find(pockets(:)== pock));
            tmpMass(fullIndex) = sum(rolledOutPermutationMatrix(pockets(:)== pock));
            fullIndex          = fullIndex + 1;
        end
        try
            threshLenUB(i)  = max(tmpLen);
            threshMassUB(i) = max(tmpMass);
        catch
            threshLenUB(i)  = 0;
            threshMassUB(i) = 0;
        end
        fullIndex = 1; tmpLen = []; tmpMass = [];
    end
    
    threshLenLB = [];
    threshMassLB = [];
    fullIndex = 1;
    
    % find living islands of below
    for i = 1:size(permutationMatrix, dim)
        
        if dim == 2
            rolledOutPermutationMatrix = squeeze(permutationMatrix(:, i)); rolledOutPermutationMatrix = rolledOutPermutationMatrix(:);
            [pockets, numPockets]      = bwlabel(belowBounds(:,i));
        elseif dim == 3
            rolledOutPermutationMatrix = squeeze(permutationMatrix(:, :, i)); rolledOutPermutationMatrix = rolledOutPermutationMatrix(:);
            [pockets, numPockets]      = bwlabel(belowBounds(:,:,i));
        end

        for pock = 1:numPockets
            tmpLen(fullIndex)  = length(find(pockets(:)== pock));
            tmpMass(fullIndex) = sum(rolledOutPermutationMatrix(pockets(:)== pock) * -1); % we multiply with -1 such that we can take the max; abs would distort the null
            fullIndex          = fullIndex + 1;
        end
        try
            threshLenLB(i)  = max(tmpLen);
            threshMassLB(i) = max(tmpMass);
        catch
            threshLenLB(i)  = 0;
            threshMassLB(i) = 0;
        end
        fullIndex = 1; tmpLen = []; tmpMass = [];
    end

    % find cluster size of LB and UB

    % these are thresholds
    clusterLenUB = prctile(threshLenUB, thresh);
    clusterLenLB = prctile(threshLenLB, thresh);
    
	clusterMassUB = prctile(threshMassUB, thresh);
    clusterMassLB = prctile(threshMassLB, thresh); % note this is inverted (16/6/23) such that we are looking at positive differences

    % these are the true datapoints
    trueUB = trueMatrix > threshUB;

    % we need to recompute threshLB here but for the inverted case
    threshLB = prctile(permutationMatrix * -1, thresh, dim);
    trueLB   = trueMatrix*-1 > threshLB;
    
    % how many of them have islands exceeding the cluster size of interest?
    [truePockets, numPockets] = bwlabel(trueUB);

    trueLenUB = [];
    trueMassUB = [];
    trueIndxUB = {};
    
    for i = 1:numPockets
    	rolledOutTruePockets = truePockets(:);
        indices = find(rolledOutTruePockets == i);

        trueLenUB(i)  = length(indices);
        trueMassUB(i) = sum(trueMatrix(indices)); % removed abs here 16/6/23 % sum(abs(trueMatrix(indices)))
        trueIndxUB{i} = indices;
    end
    
	[truePockets, numPockets] = bwlabel(trueLB);
    
	trueLenLB = [];
    trueMassLB = [];
	trueIndxLB = {};
    
    for i = 1:numPockets
    	rolledOutTruePockets = truePockets(:);
        indices = find(rolledOutTruePockets == i);
        
        trueLenLB(i)  = length(indices);
        trueMassLB(i) = sum(trueMatrix(indices) * -1); % removed abs here 16/6/23 % sum(abs(trueMatrix(indices))). note this is inverted
        trueIndxLB{i} = indices;
        
    end

    out            = struct;
    out.clustLenUB = clusterLenUB;
    out.clustLenLB = clusterLenLB;
	
    out.clustMassUB = clusterMassUB;
    out.clustMassLB = clusterMassLB; % note this is inverted
    
    out.trueUB      = trueUB;
    out.trueLB      = trueLB;
	
    out.trueLenUB   = trueLenUB;
    out.trueLenLB   = trueLenLB;
    
    out.trueMassUB = trueMassUB;
    out.trueMassLB = trueMassLB; % note this is inverted now, we ran the perm test on the inverted case to keep the same boolean thresholding logic
    
    out.trueIndxUB = trueIndxUB;
    out.trueIndxLB = trueIndxLB;
    
    out.threshUB = threshUB;
    out.threshLB = threshLB;

    indx_sig_ub = find(trueLenUB > clusterLenUB);
    indx_sig_lb = find(trueLenLB > clusterLenLB);

    indx_sig_mass_ub = find(trueMassUB > clusterMassUB);
    indx_sig_mass_lb = find(trueMassLB > clusterMassLB);

    if dim==2

        tmp_len_ub = nan(1, size(trueUB, 1));
        tmp_len_lb = nan(1, size(trueLB, 1));
    
        tmp_mass_ub = nan(1, size(trueUB, 1));
        tmp_mass_lb = nan(1, size(trueLB, 1));    

    elseif dim == 3
        tmp_len_ub = nan(size(trueUB));
        tmp_len_lb = nan(size(trueLB));
    
        tmp_mass_ub = nan(size(trueUB));
        tmp_mass_lb = nan(size(trueLB));    
    end

    try
        for e = indx_sig_ub
            tmp_len_ub(trueIndxUB{e}) = 1;
        end
    catch
    end

    try
        for e = indx_sig_lb
            tmp_len_lb(trueIndxLB{e}) = 1;
        end
    catch
    end

    try
        for e = indx_sig_mass_ub
            tmp_mass_ub(trueIndxUB{e}) = 1;
        end
    catch
    end

    try
        for e = indx_sig_mass_lb
            tmp_mass_lb(trueIndxLB{e}) = 1;
        end
    catch
    end

    out.survived_len_UB = tmp_len_ub;
    out.survived_len_LB = tmp_len_lb;

    out.survived_mass_UB = tmp_mass_ub;
    out.survived_mass_LB = tmp_mass_lb;    

end


