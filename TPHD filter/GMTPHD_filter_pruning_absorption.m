function  [weights_o,means_o,covs_o,t_ini_o,length_o,Ncom_o,logical_activas_o,means_k_old_o]=GMTPHD_filter_pruning_absorption(weights_u,means_u,covs_u,t_ini_u,length_u,logical_activas_u,T_pruning,Ncom_max,T_absorption,means_k_old,Lscan)

%Pruning and absorption algorithm for TPHD filter
%Author: Angel F. Garcia-Fernandez

weights_u(~logical_activas_u)=0;
logical_candidates=weights_u>T_pruning;

weights_o=zeros(size(weights_u));
means_o=zeros(size(means_u));
covs_o=zeros(size(covs_u));
means_k_old_o=zeros(size(means_k_old));
index_o=1;
t_ini_o=zeros(size(t_ini_u));
length_o=zeros(size(length_u));


while(sum(logical_candidates)>0)
    [~,index_max]=max(weights_u(logical_candidates));
    
    index_candidates=find(logical_candidates==1);
    
    %To select the right indices for the mean
    length_u_max=length_u(index_candidates(index_max));  
    if(length_u_max<Lscan)
        index_mean_max=4*length_u_max-3:4*length_u_max;
    else
        index_mean_max=4*Lscan-3:4*Lscan;       
    end   
    mean_max=means_u(index_mean_max,index_candidates(index_max));
    
    weight_output=0;
    logical_candidates_output=logical_candidates;
    
    for i=1:length(index_candidates)
        
        %To select the right indices for the mean and covariance
        length_k_i=length_u(index_candidates(i));
        if(length_k_i<Lscan)
            index_actual_i=4*length_k_i-3:4*length_k_i;
        else
            index_actual_i=4*Lscan-3:4*Lscan;
        end
               
        mean_i=means_u(index_actual_i,index_candidates(i));
        cov_i=covs_u(index_actual_i,index_actual_i,index_candidates(i));        
        resta=mean_i-mean_max;
        mahalanobis=resta'/cov_i*resta;  
        if(mahalanobis<=T_absorption)
            logical_candidates_output(index_candidates(i))=0;
            weight_i=weights_u(index_candidates(i));
            weight_output=weight_output+weight_i;
            
        end
    end
      
    weights_o(index_o)=weight_output;
    means_o(:,index_o)=means_u(:,index_candidates(index_max));
    covs_o(:,:,index_o)=covs_u(:,:,index_candidates(index_max)); 
    means_k_old_o(:,index_o)=means_k_old(:,index_candidates(index_max));
    t_ini_o(index_o)=t_ini_u(index_candidates(index_max));
    length_o(index_o)=length_u(index_candidates(index_max));
    
    if(index_o==Ncom_max)
        break;
    end
    
    index_o=index_o+1;
    logical_candidates=logical_candidates_output;
    
end
Ncom_o=index_o-1;
logical_activas_o=false(size(logical_activas_u));
logical_activas_o(1:Ncom_o)=1;


