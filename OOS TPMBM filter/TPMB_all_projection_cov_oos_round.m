function filter_upd_pmb=TPMB_all_projection_cov_oos_round(filter_upd,T_alive,Nx,k,Lscan)

%This function obtains the best fitting trajectory PMB to the updated
%trajectory PMBM based on KLD minimisation

%This version includes the covariance matrix for hypotheses that consider dead trajectories

%It is also tailored for the projection at OOS time. The version without
%oos is tailored to a TPMB filtering recursion with round OOS update

%Author: Angel F. Garcia-Fernandez


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
            
            
            
            flag_Bernoulli=0; %Flag=1 if it is already in Bernolli form
            if(length(filter_upd.tracks{i}.eB)==1)
                %We do not need to do anything as there is only one Bernoulli
                flag_Bernoulli=1;
            end
            
            if(flag_Bernoulli==0)
                %We do the TPMB projection
                %t_b=filter_upd.tracks{i}.t_b;
                t_ini=filter_upd.tracks{i}.t_ini;
                
                %We go through all global hypotheses and update the
                %probability of existence and the last state
                              
                
                length_m=length(filter_upd.tracks{i}.meanB{1});
                
                mean_fused=zeros(length_m,1);
                cov_fused=zeros(min(length_m,Nx*Lscan));
                prob_length_fused=zeros(size(filter_upd.tracks{i}.prob_length{1}));
                
                mean_past_fused=cell(k-t_ini,1);
                cov_past_fused=cell(Lscan-2,1);
                r_fused_past_k=zeros(Lscan-2,1);
                
                
                %Size of mean_i might not match size of cov_i due to the Lscan
                %approximation, so we obtain the right dimensions
                l_vector_i=size(cov_fused,1);
                
                r_fused=0;
                r_fused_k=0;
                
                Nhyp_i=length(filter_upd.tracks{i}.eB);
                for j=1:Nhyp_i
                    indeces_hyp_j=globHyp(:,i)==j;
                    weight_j=sum(weights(indeces_hyp_j));
                    if(weight_j>0)
                        prob_length_j=filter_upd.tracks{i}.prob_length{j};
                        prob_alive_j=prob_length_j(1);
                        
                        prod_we=weight_j*filter_upd.tracks{i}.eB(j);
                        r_fused=r_fused+prod_we;
                        r_fused_k=r_fused_k+prod_we*prob_alive_j;
                        
                        cov_i=filter_upd.tracks{i}.covB{j};
                        
                        mean_i=filter_upd.tracks{i}.meanB{j};
                        
                    
                        
                        mean_fused=mean_fused+prod_we*prob_alive_j*mean_i;
                        mean_i_Lscan=mean_i(end-l_vector_i+1:end);
                        cov_fused=cov_fused+prod_we*prob_alive_j*(cov_i+(mean_i_Lscan*mean_i_Lscan'));
                        
                        
                        
                        %The following indexing of prob_length_fused is the difference of the round version
                        prob_length_fused(1:length(prob_length_j))=prob_length_fused(1:length(prob_length_j))+prod_we*prob_length_j; 
                        
                        
                        for p=1:length(filter_upd.tracks{i}.mean_past{j})
                            
                            if(p<Lscan-1)
                                weight_jp=prod_we*prob_length_j(p+1);
                                r_fused_past_k(p)=r_fused_past_k(p)+weight_jp;
                                if(isempty(mean_past_fused{p}))
                                    %We initialise mean_past_fused{p} and
                                    %cov_past_fused{p}
                                    
                                    mean_past_fused{p}=weight_jp*filter_upd.tracks{i}.mean_past{j}{p};
                                    
                                    size_cov=size(filter_upd.tracks{i}.cov_past{j}{p},1);
                                    mean_jp_Lscan=filter_upd.tracks{i}.mean_past{j}{p}(end-size_cov+1:end);
                                    
                                    cov_past_fused{p}=weight_jp*(filter_upd.tracks{i}.cov_past{j}{p}+(mean_jp_Lscan*mean_jp_Lscan'));
                                    
                                else
                                    %We sum the new entries
                                    mean_past_fused{p}=mean_past_fused{p}+weight_jp*filter_upd.tracks{i}.mean_past{j}{p};
                                    size_cov=size(filter_upd.tracks{i}.cov_past{j}{p},1);
                                    mean_jp_Lscan=filter_upd.tracks{i}.mean_past{j}{p}(end-size_cov+1:end);
                                    
                                    cov_past_fused{p}=cov_past_fused{p}+weight_jp*(filter_upd.tracks{i}.cov_past{j}{p}+(mean_jp_Lscan*mean_jp_Lscan'));
                                    
                                end
                                
                                
                            else
                                
                                %The past trajectories beyond Lscan window
                                %are all the same.
                                mean_past_fused{p}=filter_upd.tracks{i}.mean_past{j}{p};
                                
                            end
                            
                            
                        end
                        
                        
                        
                    end
                end
                mean_fused=mean_fused/r_fused_k;   %When we consider all trajectories, we use r_fused_k, not r_fused
                mean_fused_Lscan=mean_fused(end-l_vector_i+1:end);
                cov_fused=cov_fused/r_fused_k-(mean_fused_Lscan*mean_fused_Lscan');
                prob_length_fused=prob_length_fused/r_fused; %Where we have r_fused
                
                
                for p=1:Lscan-2
                    if(r_fused_past_k(p)>0)
                        mean_past_fused{p}=mean_past_fused{p}/r_fused_past_k(p);
                        size_cov=size(cov_past_fused{p},1);
                        mean_past_fused_p_Lscan=mean_past_fused{p}(end-size_cov+1:end);
                        cov_past_fused{p}=cov_past_fused{p}/r_fused_past_k(p)-(mean_past_fused_p_Lscan*mean_past_fused_p_Lscan');
                    end
                end
                
                
                
                filter_upd_pmb.tracks{index_output}.eB=r_fused;
                filter_upd_pmb.tracks{index_output}.meanB{1}=mean_fused;
                filter_upd_pmb.tracks{index_output}.covB{1}=cov_fused;
                filter_upd_pmb.tracks{index_output}.prob_length{1}=prob_length_fused;
                
                %If the trajectory is alive with probability one (due to
                %output of Murties algorithm, there is not mean_past
                if(sum(prob_length_fused(2:end,:))>0)
                    filter_upd_pmb.tracks{index_output}.mean_past{1}=mean_past_fused;
                    filter_upd_pmb.tracks{index_output}.cov_past{1}=cov_past_fused;
                    
                else
                    filter_upd_pmb.tracks{index_output}.mean_past{1}=cell(0,1);
                    filter_upd_pmb.tracks{index_output}.cov_past{1}=cell(0,1);
                    
                end
                
                
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
                %might be empty)
                indeces=globHyp(:,i)>0;
                weight_i=sum(weights(indeces));
                
                filter_upd_pmb.tracks{index_output}.eB= weight_i*filter_upd_pmb.tracks{index_output}.eB;
                
            end
            
        end
    end
    
    %The PMB output only has one global hypothesis
    filter_upd_pmb.globHypWeight=1;
    filter_upd_pmb.globHyp=ones(1,index_output);
    
else
    %There is only one global hypotheses. It is already in PMB form.
    filter_upd_pmb=filter_upd;
    
end
