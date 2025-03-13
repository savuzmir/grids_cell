function perm_histogram(perm_inp, true_inp, hist_aux, plot_aux)
% perm_inp is a vector of permuted values
% true_inp is a value or vector of true values 
% hist_aux and plot_aux are cells containing cosmetic information you would
% put in 
if nargin == 3
    histogram(perm_inp, hist_aux{:}); hold on
elseif nargin == 4
    histogram(perm_inp, hist_aux{:});  hold on
    plot(true_inp{:}, plot_aux{:});
else
    histogram(perm_inp);  hold on
    plot(true_inp{:});   
end

hold off 