function X_estimate=GMPHD_estimation2(weights_u, means_u)

%Estimator based on Section 9.5.4.4 in R. P. S. Mahler, Advances in Statistical Multisource-Multitarget Information Fusion. Artech House, 2014
%The number of detection depends on the overall weight

N_detections=round(sum(weights_u));

if(N_detections>0)
    X_estimate=zeros(size(means_u,1)*N_detections,1);
    
    
    weights_u_copy=weights_u;
    for i=1:N_detections
        [m,index]=max(weights_u_copy);
        weights_u_copy(index)=weights_u_copy(index)-1;     
        X_estimate(4*i-3:4*i)=means_u(:,index);
     
    end
    
else
    X_estimate=[];
end
