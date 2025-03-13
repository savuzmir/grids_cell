function out = run_cluster_based_perm_test(inp, test_type, nperms, perm_type, stat_crit, dims)
% takes inp [nObservations x nTimepoints x nConditions] matrix and performs a paired or unpaired ttest across num_perms
% (inp, test_type, nperms, stat_crit, dims)
% keyboard
out_tvals = compute_perm_stats(inp, test_type, perm_type, nperms);

emp_vals          = out_tvals(:, 1);     % assumes 1st np is the true empirical mat
perm_vals         = out_tvals(:, 2:end); % assumes 2:end np is theoretical mat  

out_cluster_stats = compute_cluster_stats(emp_vals, perm_vals, stat_crit, dims);

out = struct;
out.orig_tvals    = out_tvals;
out.cluster_stats = out_cluster_stats;

end