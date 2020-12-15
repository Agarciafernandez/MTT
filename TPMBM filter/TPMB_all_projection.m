function filter_upd_pmb=TPMB_all_projection(filter_upd,T_alive,Nx,k,Lscan)

%This function obtains the best fitting trajectory PMB to the updated
%trajectory PMBM based on KLD minimisation

%Author: Angel F. Garcia-Fernandez
%Copyright (c) 2020, Angel F. Garcia-Fernandez


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
            alive_trajectory=0;
            
            if(length(filter_upd.tracks{i}.eB)==1)
                %We do not need to do anything as there is only one Bernoulli
                flag_Bernoulli=1;
            else
                
                %We check if the Bernoulli component is alive in any of the
                %hypotheses
                for j=1:length(weights)
                    index_hyp_i_j=globHyp(j,i); %Track index
                    if(index_hyp_i_j>0)
                        prob_length_j=filter_upd.tracks{i}.prob_length{index_hyp_i_j};
                        prob_alive_j=prob_length_j(1);
                        if(prob_alive_j>T_alive)
                            alive_trajectory=1;
                            break;
                        end
                    end
                end
            end
            
            if(flag_Bernoulli==0 && alive_trajectory)
                %We need to do the projection
                t_b=filter_upd.tracks{i}.t_b;
                t_ini=filter_upd.tracks{i}.t_ini;
                prob_length_fused=zeros(k-t_ini+1,1);
                
                %We go through all global hypotheses and update the
                %probability of existence and the last state
                
                mean_fused=zeros(Nx*(k-t_b+1),1);
                cov_fused=zeros(Nx*min(k-t_b+1,Lscan));
                
                
                %Size of mean_i might not match size of cov_i due to the Lscan
                %approximation, so we obtain the right dimensions
                l_vector_i=size(cov_fused,1);
                
                r_fused=0;
                r_fused_k=0;
                
                Nhyp_i=length(filter_upd.tracks{i}.eB);
                for p=1:Nhyp_i
                    indeces_hyp_p=globHyp(:,i)==p;
                    weight_p=sum(weights(indeces_hyp_p));
                    if(weight_p>0)
                        prob_length_j=filter_upd.tracks{i}.prob_length{p};
                        prob_alive_j=prob_length_j(1);
                        
                        prod_we=weight_p*filter_upd.tracks{i}.eB(p);
                        r_fused=r_fused+prod_we;
                        r_fused_k=r_fused_k+prod_we*prob_alive_j;
                        
                        cov_i=filter_upd.tracks{i}.covB{p};
                        
                        mean_i=filter_upd.tracks{i}.meanB{p};
                        
                        mean_fused=mean_fused+prod_we*prob_alive_j*mean_i;
                        mean_i_Lscan=mean_i(end-l_vector_i+1:end);
                        cov_fused=cov_fused+prod_we*prob_alive_j*(cov_i+(mean_i_Lscan*mean_i_Lscan'));
                        
                        if(prob_alive_j<1)
                            %All the past trajectories are the same in hypotheses in which prob_alive_j<1
                            mean_past_fused=filter_upd.tracks{i}.mean_past{p};
                            %We update prob_length_fused
                            prob_length_fused=prob_length_fused+prod_we*prob_length_j;
                            
                        else
                            prob_length_fused(1)=prob_length_fused(1)+prod_we*prob_alive_j;
                        end
                        
                        
                    end
                end
                mean_fused=mean_fused/r_fused_k;   %When we consider all trajectories, we use r_fused_k, not r_fused
                mean_fused_Lscan=mean_fused(end-l_vector_i+1:end);
                cov_fused=cov_fused/r_fused_k-(mean_fused_Lscan*mean_fused_Lscan');
                prob_length_fused=prob_length_fused/r_fused; %Where we have r_fused
                
                
                
                
                
                filter_upd_pmb.tracks{index_output}.eB=r_fused;
                filter_upd_pmb.tracks{index_output}.meanB{1}=mean_fused;
                filter_upd_pmb.tracks{index_output}.covB{1}=cov_fused;
                filter_upd_pmb.tracks{index_output}.prob_length{1}=prob_length_fused;
                
                %If the trajectory is alive with probability one (due to
                %output of Murties algorithm, there is not mean_past
                if(sum(prob_length_fused(2:end,:))>0)
                    filter_upd_pmb.tracks{index_output}.mean_past{1}=mean_past_fused;
                else
                    filter_upd_pmb.tracks{index_output}.mean_past{1}=cell(0,1);
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
                %might be empty
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
