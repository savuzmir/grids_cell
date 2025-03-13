function [unique_fitted_hex, unique_orig_hex, unique_hex_diff, rel_sess] = compute_session_level_stats(animal, session, fitted, original,brain_region,rel_rg, nTR)

indx = 1;

    for an = 0:1
    
        animal_inds = find(animal==an);     
        animal_sess = session(animal_inds);
        
        fitted_hex = fitted(:, animal_inds, :);
        orig_hex   = original(:, animal_inds, :);
        hex_diff   = fitted_hex - orig_hex;
        
        animal_region = brain_region(animal_inds);
        session_IDs = unique(animal_sess);
    
        for sess = 1:length(session_IDs)
            
            tmpBuffer= [];
            
            curr_sess = session_IDs(sess);
           
            for rg = rel_rg
                tmpBuffer = [tmpBuffer; find((animal_sess == curr_sess) &  (animal_region == rg))'];
            end
    
            fitted_hex_sess = [];
            orig_hex_sess   = [];
            hex_diff_sess    = [];
            
            % this gives us session averages for our selected regions in rg
            fitted_hex_sess(:, :, sess) = squeeze(nanmean(fitted_hex(:, tmpBuffer', :), 2));
            orig_hex_sess(:, :, sess)   = squeeze(nanmean(orig_hex(:, tmpBuffer', :), 2));
            hex_diff_sess(:, :, sess)   = squeeze(nanmean(hex_diff(:, tmpBuffer', :), 2));
      
            % this gives session means
            unique_fitted_hex(:, :, indx) = fitted_hex_sess(:, :, sess);
            unique_orig_hex(:, :, indx) = orig_hex_sess(:, :, sess);
            unique_hex_diff(:, :, indx) = hex_diff_sess(:, :, sess);
    
            rel_sess(indx) = nanmean(nTR(tmpBuffer));
            
            indx = indx + 1;
        end  
    end

end