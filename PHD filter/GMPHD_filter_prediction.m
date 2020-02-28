function   [weights_k, means_k,covs_k,Ncom_k,logical_active_k]=GMPHD_filter_prediction(weights_u, means_u,covs_u,Ncom_u,logical_active_u,F,Q,p_s)

%Prediction for the GMPHD filter

Ncom_k=Ncom_u;
logical_active_k=logical_active_u;

index_active_u=find(logical_active_u==1);

weights_k=zeros(size(weights_u));
means_k=zeros(size(means_u));
covs_k=zeros(size(covs_u));

for i=1:Ncom_u
    mean_i=means_u(:,index_active_u(i));
    cov_i=covs_u(:,:,index_active_u(i));
    weight_i=weights_u(index_active_u(i));
    
    mean_i_pred=F*mean_i;
    cov_i_pred=F*cov_i*F'+Q;
    weight_i_pred=weight_i*p_s;
    
    weights_k(index_active_u(i))=weight_i_pred;
    means_k(:,index_active_u(i))=mean_i_pred;
    covs_k(:,:,index_active_u(i))=cov_i_pred;
   
end
