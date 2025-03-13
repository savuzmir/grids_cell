Code for reproducing figures in 'A cognitive map for value-guided choice in ventromedial prefrontal cortex'. 
The code is structured such that figure panels from individual figures can be regenerated from the 'main_grid' function:

Figure 1:
figure1_BehaviourValue();
figure1_supplements_warping();

Figure 2:
figure2_lfpGrid();
figure2_lfpGrid_supplements();

Figure 3: 
figure3_gridRealignment();

Figure 4:
figure4_GridNeurons();
figure4_GridNeurons_supplements(); 
ofc_post_revision_supplements();

Figure 5:
figure5_ripples();
figure5_rippleSupplements(); 
nnmf_example_supplements(); 

Note that loading in data for some of the functions/figures consumes a lot of RAM so may cause loading issues.
