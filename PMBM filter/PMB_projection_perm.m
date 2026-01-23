function filter_upd_pmb=PMB_projection_perm(filter_upd,Perm_hyp)

%PMB projection using a permutation of the auxiliary variables for each
%global hypothesis
%This is part of the variational PMB filter (VPMB_projection.m)
%Author: Ángel García-Fernández

%Weights of global hypotheses
weights=filter_upd.globHypWeight;

Nx=size(filter_upd.meanPois,1);

Nglob_hyp=length(weights);

if(Nglob_hyp>1)
    %The PPP does not change
    filter_upd_pmb.weightPois=filter_upd.weightPois;
    filter_upd_pmb.meanPois=filter_upd.meanPois;
    filter_upd_pmb.covPois=filter_upd.covPois;

    Ntracks=length(filter_upd.tracks);
    globHyp=filter_upd.globHyp;


    %We go through each Bernoulli component

    for i=1:Ntracks

        %Fused parameters of the Bernoulli
        r_fused=0;
        mean_fused=zeros(Nx,1);
        cov_fused=zeros(Nx,Nx);

        %We go through all global hypotheses
        for p=1:Nglob_hyp


         

            weight_p=weights(p);
            Perm_hyp_i_p=Perm_hyp(i,p);

            index_hyp_p=globHyp(p,Perm_hyp_i_p); %We choose the local hypotheses for the permuted Bernoulli (index Perm_hyp_i_p) in this global hypothesis

            if(index_hyp_p>0)
                %This Bernoulli exists in this local hypothesis (otherwise,
                %we do not do anything)
                eB=filter_upd.tracks{Perm_hyp_i_p}.eB(index_hyp_p);
                meanB=filter_upd.tracks{Perm_hyp_i_p}.meanB{index_hyp_p};
                covB=filter_upd.tracks{Perm_hyp_i_p}.covB{index_hyp_p};

                prod_we=weight_p*eB;
                r_fused=r_fused+prod_we;

                mean_fused=mean_fused+prod_we*meanB;
                cov_fused=cov_fused+prod_we*(covB+(meanB*meanB'));

            end

        end

        filter_upd_pmb.tracks{i}.eB=min(r_fused,1); %We use this min to improve numerical accuracy

        if(r_fused>0)
            mean_fused=mean_fused/r_fused;
            cov_fused=cov_fused/r_fused-(mean_fused*mean_fused');

            filter_upd_pmb.tracks{i}.meanB{1}=mean_fused;
            filter_upd_pmb.tracks{i}.covB{1}=cov_fused;
        else
            filter_upd_pmb.tracks{i}.meanB{1}=NaN(Nx,1);
            filter_upd_pmb.tracks{i}.covB{1}=NaN(Nx,Nx);

        end

        %The rest of the parameters do not change
        filter_upd_pmb.tracks{i}.t_ini=filter_upd.tracks{i}.t_ini;

        filter_upd_pmb.tracks{i}.aHis{1}=0; %The data association variable is meaningless in the PMB filters but it is required so that the code works


    end

    %The MB output only has one global hypothesis
    filter_upd_pmb.globHypWeight=1;
    filter_upd_pmb.globHyp=ones(1,Ntracks);

else
    %There is only one global hypotheses. It is already in PMB form.
    filter_upd_pmb=filter_upd;

end
