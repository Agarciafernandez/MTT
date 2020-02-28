function  [weights_o, means_o,covs_o,Ncom_o,logical_actives_o]=GMPHD_filter_pruning(weights_u, means_u,covs_u,logical_active_u,Ncom_max,T_pruning,T_merging)
%Pruning and merging of PHD components in GMPHD filter

weights_u(~logical_active_u)=0;
logical_candidates=weights_u>T_pruning;

weights_o=zeros(size(weights_u));
means_o=zeros(size(means_u));
covs_o=zeros(size(covs_u));
index_o=1;
Nx=size(means_u,1);


while(sum(logical_candidates)>0)
    
    [c,index_max]=max(weights_u(logical_candidates));
    
    index_candidates=find(logical_candidates==1);
    
    mean_max=means_u(:,index_candidates(index_max));
    
    weight_output=0;
    mean_output=zeros(Nx,1);
    cov_output=zeros(Nx,Nx);
    logical_candidates_output=logical_candidates;
    
    for i=1:length(index_candidates)
        mean_i=means_u(:,index_candidates(i));
        cov_i=covs_u(:,:,index_candidates(i));        
        resta=mean_i-mean_max;
        mahalanobis=resta'/cov_i*resta;
        
      
        
        if(mahalanobis<T_merging)
            logical_candidates_output(index_candidates(i))=0;
            weight_i=weights_u(index_candidates(i));
            weight_output=weight_output+weight_i;
            mean_output=mean_output+weight_i*mean_i;
            cov_output=cov_output+weight_i*(cov_i+mean_i*mean_i');
        end
    end
    
    mean_output=mean_output/weight_output;
    cov_output=cov_output/weight_output-mean_output*mean_output';
    
    weights_o(index_o)=weight_output;
    means_o(:,index_o)=mean_output;
    covs_o(:,:,index_o)=(cov_output+cov_output')/2;   
    
    if(index_o==Ncom_max)
        break;
    end
    
    index_o=index_o+1;
    logical_candidates=logical_candidates_output;
    
end
Ncom_o=index_o-1;
logical_actives_o=false(size(logical_active_u));
logical_actives_o(1:Ncom_o)=1;


