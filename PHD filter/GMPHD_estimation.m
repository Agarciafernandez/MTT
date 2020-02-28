function X_estimate=GMPHD_estimation(weights_u, means_u)

%Estimator based on B.-N. Vo and W.-K. Ma, “The Gaussian mixture probability hypothesis density filter,” IEEE Trans. Signal Process., vol. 54, no. 11, pp. 4091–4104, Nov. 2006.

detections=round(weights_u);
N_detections=sum(detections);
if(N_detections>0)
    X_estimate=zeros(size(means_u,1),N_detections);
    index_min=1;  
    while(N_detections>0)
        N_detections_i=sum(detections>0);
        X_estimate(:,index_min:N_detections_i+index_min-1)=means_u(:,detections>0);
        
        detections=detections-(detections>0);
        N_detections=sum(detections);
        index_min=N_detections_i+index_min;          
    end    
    X_estimate=X_estimate(:);
else
    X_estimate=[];
end