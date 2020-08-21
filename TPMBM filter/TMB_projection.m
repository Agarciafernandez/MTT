function filter_upd_pmb=TMB_projection(filter_upd)

%This function obtains the best fitting trajectory MB to the updated
%trajectory MBM based on KLD minimisation
%It does not consider the Poisson part, compared to the TPMB projection 

%Faster version of the projection


%Weights of global hypotheses
weights=filter_upd.globHypWeight;


if(length(weights)>1)
    
    Ntracks=length(filter_upd.tracks);
    globHyp=filter_upd.globHyp;
    
    %The PMB output only has one global hypothesis
    filter_upd_pmb.globHypWeight=1;
    filter_upd_pmb.globHyp=ones(1,Ntracks);
    
    
    %We go through each Bernoulli component
    
    for i=1:Ntracks
        
        if(length(filter_upd.tracks{i}.eB)>1)
            
            
            %We go through all global hypotheses
            r_fused=0;
            mean_fused=zeros(size(filter_upd.tracks{i}.meanB{1}));
            cov_fused=zeros(size(filter_upd.tracks{i}.covB{1}));
            
            %Size of mean_i might not match size of cov_i due to the Lscan
            %approximation, so we obtain the right dimensions
            l_vector_i=size(cov_fused,1);
            
            
            %             Option 1 for mean and covariances: moment matching
            
            
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
                    
                   % cov_fused=cov_fused+prod_we*cov_i;
                    
                end
            end
            mean_fused=mean_fused/r_fused;
            mean_fused_Lscan=mean_fused(end-l_vector_i+1:end);
            cov_fused=cov_fused/r_fused-(mean_fused_Lscan*mean_fused_Lscan');
            
            
            cov_fused=(cov_fused+cov_fused')/2;
            
%                         mean_fused=mean_fused/r_fused;
%                         mean_fused_Lscan=mean_fused(end-l_vector_i+1:end);
%                          cov_fused=cov_fused/r_fused;
            
            
            %Option 2 for mean and covariances: compute the existence as
            %before but select component with highest weight (according to
            %prod_we)
            
            %             Nhyp_i=length(filter_upd.tracks{i}.eB);
            %             prod_we_max=0;
            %
            %             for p=1:Nhyp_i
            %                 indeces_hyp_p=globHyp(:,i)==p;
            %                 weight_p=sum(weights(indeces_hyp_p));
            %                 if(weight_p>0)
            %                     prod_we=weight_p*filter_upd.tracks{i}.eB(p);
            %                     r_fused=r_fused+prod_we;
            %                     if(prod_we>prod_we_max)
            %                         mean_i=filter_upd.tracks{i}.meanB{p};
            %                         cov_i=filter_upd.tracks{i}.covB{p};
            %                         mean_fused=mean_i;
            %                         cov_fused=cov_i;
            %                         prod_we_max=prod_we;
            %                     end
            %                 end
            %             end
            
            
            
            filter_upd_pmb.tracks{i}.eB=r_fused;
            filter_upd_pmb.tracks{i}.meanB{1}=mean_fused;
            filter_upd_pmb.tracks{i}.covB{1}=cov_fused;
            
            %The rest of the parameters do not change
            filter_upd_pmb.tracks{i}.t_ini=filter_upd.tracks{i}.t_ini;
            filter_upd_pmb.tracks{i}.t_b=filter_upd.tracks{i}.t_b;
            filter_upd_pmb.tracks{i}.length=filter_upd.tracks{i}.length;
            
            %The data association variable is meaningless in the TPMB filter but is required so that
            %the code works. We nevertheless take filter_upd.tracks{i}.aHis{i}(1), as it captures the "label" together with initial time step
            
            filter_upd_pmb.tracks{i}.aHis{1}=filter_upd.tracks{i}.aHis{1}(1);
            
        else
            %The Bernoulli component is already Bernoulli
            filter_upd_pmb.tracks{i}=filter_upd.tracks{i};
            
        end
        
    end
    
else
    %There is only one global hypotheses. It is already in PMB form.
    filter_upd_pmb=filter_upd;
    
end
