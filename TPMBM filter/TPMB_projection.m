function filter_upd_pmb=TPMB_projection(filter_upd)

%This function obtains the best fitting trajectory PMB to the updated
%trajectory PMBM based on KLD minimisation with auxiliary variables


%Weights of global hypotheses
weights=filter_upd.globHypWeight;


if(length(weights)>1)
    %The PPP does not change
    filter_upd_pmb.Pois=filter_upd.Pois;
    
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
                
                %Size of mean_i might not match size of cov_i due to the Lscan
                %approximation, so we obtain the right dimensions
                l_vector_i=size(cov_fused,1);
                
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
                        
                        mean_i_Lscan=mean_i(end-l_vector_i+1:end);
                        
                        cov_fused=cov_fused+prod_we*(cov_i+(mean_i_Lscan*mean_i_Lscan'));
                        
                        
                        
                    end
                    
                end
                mean_fused=mean_fused/r_fused;
                mean_fused_Lscan=mean_fused(end-l_vector_i+1:end);
                cov_fused=cov_fused/r_fused-(mean_fused_Lscan*mean_fused_Lscan');
                
                
                filter_upd_pmb.tracks{index_output}.eB=r_fused;
                filter_upd_pmb.tracks{index_output}.meanB{1}=mean_fused;
                filter_upd_pmb.tracks{index_output}.covB{1}=cov_fused;
                
                %The rest of the parameters do not change
                filter_upd_pmb.tracks{index_output}.t_ini=filter_upd.tracks{i}.t_ini;
                filter_upd_pmb.tracks{index_output}.t_b=filter_upd.tracks{i}.t_b;
                filter_upd_pmb.tracks{index_output}.length=filter_upd.tracks{i}.length;
                
                
                filter_upd_pmb.tracks{index_output}.aHis{1}=0; %The data association variable is meaningless in the TPMB filter but is required so that the code works
                
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
