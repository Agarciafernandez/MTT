function filter_upd_pmb=PMB_projection(filter_upd)

%This function obtains the best fitting PMB to the updated
% PMBM based on KLD minimisation with auxiliary variable [1].
% This projection gives rise to the track-oriented PMB filter in [2].

%[1] Á. F. García-Fernández, L. Svensson, J. L. Williams, Y. Xia and K. Granström, 
% "Trajectory Poisson Multi-Bernoulli Filters," in IEEE Transactions on Signal Processing, vol. 68, pp. 4933-4945, 2020

%[2] J. L. Williams, "Marginal multi-bernoulli filters: RFS derivation of MHT, JIPDA, and association-based member," 
% in IEEE Transactions on Aerospace and Electronic Systems, vol. 51, no. 3, pp. 1664-1687, July 2015,


%Weights of global hypotheses
weights=filter_upd.globHypWeight;


if(length(weights)>1)
    %The PPP does not change
    filter_upd_pmb.weightPois=filter_upd.weightPois;
    filter_upd_pmb.meanPois=filter_upd.meanPois;
    filter_upd_pmb.covPois=filter_upd.covPois;
    
    Ntracks=length(filter_upd.tracks);
    globHyp=filter_upd.globHyp;
    
    index_output=0;
    %We go through each Bernoulli component
    
    for i=1:Ntracks
        
        %We do not do anything if this coresponds to a Bernoulli with
        %existence probaility =0 (e.g. clutter originated Bernoulli)
        if(length(filter_upd.tracks{i}.eB)==1 && filter_upd.tracks{i}.eB==0)
            
            %We do not do anything
            
        else
            index_output=index_output+1;
            
            if(length(filter_upd.tracks{i}.eB)>1)
                
                %We go through all global hypotheses
                r_fused=0;
                mean_fused=zeros(size(filter_upd.tracks{i}.meanB{1}));
                cov_fused=zeros(size(filter_upd.tracks{i}.covB{1}));
                
                Nhyp_i=length(filter_upd.tracks{i}.eB);
                for p=1:Nhyp_i
                    indeces_hyp_p=globHyp(:,i)==p;
                    weight_p=sum(weights(indeces_hyp_p));
                    if(weight_p>0)
                        
                        prod_we=weight_p*filter_upd.tracks{i}.eB(p);
                        r_fused=r_fused+prod_we;
                        cov_i=filter_upd.tracks{i}.covB{p};
                        
                        %mean_i_Lscan=filter_upd.tracks{i}.meanB{index_hyp_i_j}(end-l_vector_i+1:end);
                        
                        mean_i=filter_upd.tracks{i}.meanB{p};
                        mean_fused=mean_fused+prod_we*mean_i;
                        cov_fused=cov_fused+prod_we*(cov_i+(mean_i*mean_i'));
                    end
                end
                mean_fused=mean_fused/r_fused;
                cov_fused=cov_fused/r_fused-(mean_fused*mean_fused');
                
                
                filter_upd_pmb.tracks{index_output}.eB=r_fused;
                filter_upd_pmb.tracks{index_output}.meanB{1}=mean_fused;
                filter_upd_pmb.tracks{index_output}.covB{1}=cov_fused;
                
                %The rest of the parameters do not change
                filter_upd_pmb.tracks{index_output}.t_ini=filter_upd.tracks{i}.t_ini;
                
                %We keep the data association variable at the time step
                %when the i-th Bernoulli component was created. This auxiliary variable
                %can be used for sequential track formation, see for
                %instance PMBtarget_filter_tracks_all.m
                filter_upd_pmb.tracks{index_output}.aHis{1}= filter_upd.tracks{i}.aHis{1}(1); 
                
            else
                
                %The Bernoulli component is already Bernoulli
                filter_upd_pmb.tracks{index_output}=filter_upd.tracks{i};
                
                %In the PMB form, we must recalculate the existence
                %probability according to global weight (as some hypotheses
                %might be empty
                indeces=globHyp(:,i)>0;
                weight_i=sum(weights(indeces));
                
                filter_upd_pmb.tracks{index_output}.eB= weight_i*filter_upd_pmb.tracks{index_output}.eB;
            end
        end
        
    end
    
    %The MB output only has one global hypothesis
    filter_upd_pmb.globHypWeight=1;
    filter_upd_pmb.globHyp=ones(1,index_output);
    
else
    %There is only one global hypotheses. It is already in PMB form.
    filter_upd_pmb=filter_upd;
    
end
