function pmbm_fusion=NormProductPMBs_assigned(pmb_1,pmb_2,pmbm_fusion,gating_threshold,Nprev_tracks,Nnew_tracks)

%This function takes an pmb_1 and updates its assigned Bernoulli components
%with the Bernoulli of pmb_2
%The output is added to the fused PMBM which is stored in pmbm_fusion


%Author: Ángel García-Fernández


for i=1:Nprev_tracks
    %As we are updating a Bernoulli, there is only one previous hypothesis

    pmbm_fusion.tracks{i}.t_ini=0; %This metadata is set to zero, as it loses its meaning once we do fusion

    mean_i=pmb_1.tracks{i}.meanB{1};
    cov_i=pmb_1.tracks{i}.covB{1};
    eB_i=pmb_1.tracks{i}.eB(1);

    %We go though all Bernoulli components of PMB2 ("measurements")
    for m=1:Nnew_tracks
        index_hyp=1+m;
        %Subindex m refers to the second PMB, which is acting as a "measurement"
        mean_m=pmb_2.tracks{m}.meanB{1};
        cov_m=pmb_2.tracks{m}.covB{1};
        eB_m=pmb_2.tracks{m}.eB(1);

        %Squared Mahalanobis distance for gating (see paper)
        maha=(mean_m-mean_i)'/(cov_m+cov_i)*(mean_m-mean_i);

        if(maha<gating_threshold)
            %We create the component

            [mean_prod,cov_prod,alpha_prod]=ProductTwoGaussians(mean_i,cov_i,mean_m,cov_m);

            rho=(1-eB_i)*(1-eB_m)+eB_i*eB_m*alpha_prod;

            eB_fused=eB_i*eB_m*alpha_prod/rho;

            pmbm_fusion.tracks{i}.meanB{index_hyp}=mean_prod;
            pmbm_fusion.tracks{i}.covB{index_hyp}=cov_prod;

            pmbm_fusion.tracks{i}.eB(index_hyp)=eB_fused;
            pmbm_fusion.tracks{i}.aHis{index_hyp}=0; %This metadata is set to zero

            %Update of the weights
            pmbm_fusion.tracks{i}.weightBLog_k(index_hyp)=log(rho);


        else
            pmbm_fusion.tracks{i}.weightBLog_k(index_hyp)=-Inf;
        end

    end
end






