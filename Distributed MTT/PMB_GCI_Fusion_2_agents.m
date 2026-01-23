function filter_agents_fused=PMB_GCI_Fusion_2_agents(filter_agents,omega,pmb_proj_after_f,gating_threshold_fusion,Nhyp_max_norm_prod,T_pruning,T_pruningPois,Nhyp_max,existence_threshold,T_mergingPois,NcomPois_max,max_iter_vpmb,threshold_vpmb)

%Author: Ángel García-Fernández

%This function performs the  GCI fusion rule for PMB densities for two agents explained in 

%Á. F. García-Fernández and G. Battistelli, "Distributed Poisson Multi-Bernoulli Filtering via Generalized Covariance Intersection," in IEEE Transactions on Signal Processing, vol. 74, pp. 246-257, 2026, doi: 10.1109/TSP.2026.3651805. 




%We calculate the power of two PMB, oen with omega and the second with 1-omega
power_pmb_1=CalculatePowerPMB(filter_agents{1}.filter_upd,omega);
power_pmb_2=CalculatePowerPMB(filter_agents{2}.filter_upd,1-omega);


%We calculate the normalised product of two PMBs
pmbm_fusion=NormProductPMBs(power_pmb_1,power_pmb_2,gating_threshold_fusion,Nhyp_max_norm_prod);

if(pmb_proj_after_f==1)
    %We do the TO-PMB projection
    pmb_fusion=PMB_projection(pmbm_fusion);

    %We write the fused PMB in each filter update
    filter_agents_fused{1}.filter_upd=pmb_fusion;
    filter_agents_fused{2}.filter_upd=pmb_fusion;
elseif(pmb_proj_after_f==2)
    %We do VPMB projection
    pmbm_fusion_pruned=PMBMtarget_pruning_merging(pmbm_fusion, T_pruning,T_pruningPois,Nhyp_max,existence_threshold,T_mergingPois,NcomPois_max);
    pmb_fusion=VPMB_projection(pmbm_fusion_pruned,max_iter_vpmb,threshold_vpmb);

    %We do pruning to clean the PMB (some Bernoulli components may have
    %probability of existence=0 after the VPMB projection, and these must be removed before the next update)
    pmb_fusion=PMBMtarget_pruning_merging(pmb_fusion, T_pruning,T_pruningPois,Nhyp_max,existence_threshold,T_mergingPois,NcomPois_max);

    filter_agents_fused{1}.filter_upd=pmb_fusion;
    filter_agents_fused{2}.filter_upd=pmb_fusion;


else

    filter_agents_fused{1}.filter_upd=pmbm_fusion;
    filter_agents_fused{2}.filter_upd=pmbm_fusion;

end



