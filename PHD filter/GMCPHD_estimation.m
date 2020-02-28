function X_estimate=GMCPHD_estimation(weights_u, means_u,cardinality_u)

[max_card,index_card]=max(cardinality_u);
N_detections=index_card-1;
if(N_detections>0)
    X_estimate=zeros(size(means_u,1),N_detections);    
    for i=1:N_detections
        [max_weight,index_max]=max(weights_u);
        X_estimate(:,i)=means_u(:,index_max);
        weights_u(index_max)=0;   
    end
    X_estimate=X_estimate(:);
else
    X_estimate=[];
end