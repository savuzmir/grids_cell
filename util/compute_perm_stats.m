function tvalue_stats = compute_perm_stats(inp, test_type, perm_type, num_perms)
% takes inp [nObservations x nTimepoints x nConditions] matrix and performs a paired or unpaired ttest across num_perms

    tmp_a = inp(:, :, 1);
    tmp_b = inp(:, :, 2);
    
    num_obs_a = sum(~isnan(tmp_a(:, 1)));
    num_obs_b = sum(~isnan(tmp_b(:, 1)));
    
    tmp_a = tmp_a(find(~isnan(tmp_a(:, 1))), :);
    tmp_b = tmp_b(find(~isnan(tmp_b(:, 1))), :);
    
    for np = 1:num_perms
    
       if np == 1
           if test_type == 1     % paired t-test for empirical data
               [~,~,~,tst_out] = ttest(tmp_a, tmp_b);
           elseif test_type == 2 % independent t-test for empirical data
               [~,~,~,tst_out] = ttest2(tmp_a, tmp_b);
           end
       else

           if perm_type == 1 % permutation of labels across conditions 
               tmp_a_b_stretch = [tmp_a; tmp_b];
        
               tmp_a_b_stretch = tmp_a_b_stretch(shuffle(1:(num_obs_a+num_obs_b)), :);
        
               tmp_a_b_stretch_1 = tmp_a_b_stretch(1:num_obs_a, :);
               tmp_a_b_stretch_2 = tmp_a_b_stretch(num_obs_a+1:(num_obs_a+num_obs_b), :);

           elseif perm_type == 2 % paired permutation 

               for kk = 1:size(tmp_a, 1) % assume same size - this should be fine as we would only do paired permutation for paired cases
                   if rand < .5
                       tmp_a_b_stretch_1(kk, :) = tmp_a(kk, :);
                       tmp_a_b_stretch_2(kk, :) = tmp_b(kk, :);
                   else
                       tmp_a_b_stretch_1(kk, :) = tmp_b(kk, :);
                       tmp_a_b_stretch_2(kk, :) = tmp_a(kk, :);
                   end
               end
           end

           if test_type == 1
               [~,~,~,tst_out] = ttest(tmp_a_b_stretch_1, tmp_a_b_stretch_2);
           else
               [~,~,~,tst_out] = ttest2(tmp_a_b_stretch_1, tmp_a_b_stretch_2);
           end
       end

       tvalue_stats(:, np) = tst_out.tstat;

    end

end