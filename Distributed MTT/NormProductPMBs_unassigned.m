function pmbm_fusion=NormProductPMBs_unassigned(pmb_1,pmb_2,pmbm_fusion,gating_threshold,Nprev_tracks,Nnew_tracks,Nx)

%This function takes an pmb_1 and updates its unassigned Bernoulli components with
%the PPP component of pmb_2. 
%The output is added to the fused PMBM which is stored in pmbm_fusion
%Nprev_tracks is used as an offset to store the newly created Bernoulli
%components. When we update pmb_1, this is set to Nprev_tracks, and when we
%update pmb_2 (calling the function pmb_1=pmb_2 and pmb_2=pmb_1, it should
%be set to zero.)

%Author: Ángel García-Fernández

for m=1:Nnew_tracks

    %Subindex m refers to the second PMB, which is acting as a "measurement"
    mean_m=pmb_2.tracks{m}.meanB{1};
    cov_m=pmb_2.tracks{m}.covB{1};
    eB_m=pmb_2.tracks{m}.eB(1);

    %We go through all Poisson components of PMB1 and perform gating

    indices_new_tracks=[];
    for i=1:length(pmb_1.weightPois)

        mean_i=pmb_1.meanPois(:,i);
        cov_i=pmb_1.covPois(:,:,i);


        %Squared Mahalanobis distance for gating (see paper)
        maha=(mean_m-mean_i)'/(cov_m+cov_i)*(mean_m-mean_i);

        if(maha<gating_threshold)
            %A track is created
            indices_new_tracks=[indices_new_tracks,i];

        end

    end

    if(~isempty(indices_new_tracks))
        %A new track is created for this measurement
        meanB=zeros(Nx,1);
        covB=zeros(Nx);
        weightB=0;

        for i=1:length(indices_new_tracks)


            index=indices_new_tracks(i);
            mean_i=pmb_1.meanPois(:,index);
            cov_i=pmb_1.covPois(:,:,index);

            [mean_prod,cov_prod,alpha_prod]=ProductTwoGaussians(mean_i,cov_i,mean_m,cov_m);

            weight_i=alpha_prod*pmb_1.weightPois(index);
            weightB=weightB+weight_i;

            %Gaussian mixture reduction
            meanB=meanB+weight_i*mean_prod;
            covB=covB+weight_i*cov_prod+weight_i*(mean_prod*mean_prod');
        end

        meanB=meanB/weightB;
        covB=covB/weightB-(meanB*meanB');
        rho=1-eB_m+eB_m*weightB;

        eB=eB_m*weightB/rho;



        pmbm_fusion.tracks{Nprev_tracks+m}.meanB{1}=meanB;
        pmbm_fusion.tracks{Nprev_tracks+m}.covB{1}=covB;
        pmbm_fusion.tracks{Nprev_tracks+m}.eB=eB;
        pmbm_fusion.tracks{Nprev_tracks+m}.t_ini=0;  %This metadata is set to zero, as it loses its meaning once we do fusion
        pmbm_fusion.tracks{Nprev_tracks+m}.aHis{1}=m; %This metadata is set to m (this information is used for obtaining the k-best global hypotheses, as required by previous code)
        pmbm_fusion.tracks{Nprev_tracks+m}.weightBLog_k=log(rho);

    else


        %The Bernoulli component has existence probability zero (it cannot be associated to PPP). It will be removed by pruning
        rho=1-eB_m;

        pmbm_fusion.tracks{Nprev_tracks+m}.eB=0;
        pmbm_fusion.tracks{Nprev_tracks+m}.weightBLog_k=log(rho+eps); %We add eps for numerical robustness in case eB_m=1

        pmbm_fusion.tracks{Nprev_tracks+m}.meanB{1}=zeros(Nx,1);
        pmbm_fusion.tracks{Nprev_tracks+m}.covB{1}=zeros(Nx,Nx);
        pmbm_fusion.tracks{Nprev_tracks+m}.t_ini=0;  %This metadata is set to zero, as it loses its meaning once we do fusion
        pmbm_fusion.tracks{Nprev_tracks+m}.aHis{1}=m; %This metadata is set to m (this information is used for obtaining the k-best global hypotheses, as required by previous code)

    end

end
