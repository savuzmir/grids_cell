function plot_sigline_perm(inp, perm_inp, perm_base, perm_tp, smooth_val)

plotmse(inp, [smooth_val(1), smooth_val(2)]); hold on

if perm_tp == 1
    sig_len    = find(~isnan(perm_inp.survived_len_UB));
    sig_len_lb = find(~isnan(perm_inp.survived_len_LB));
elseif perm_tp == 2
    sig_len    = find(~isnan(perm_inp.survived_mass_UB));
    sig_len_lb = find(~isnan(perm_inp.survived_mass_LB));
end

tmp_ub          = nan(1, size(inp, 2));
tmp_ub(sig_len) = perm_base + perm_base*0.05;

try
    plot(tmp_ub, 'linewidth', 4, 'color', 'k')
catch
    ...
end

tmp_lb             = nan(1, size(inp, 2));
tmp_lb(sig_len_lb) = perm_base + perm_base*0.05;

try
    plot(tmp_lb,  'linewidth', 4, 'color', 'k')
catch
    ...
end

end