function    [X_estimate,t_b_estimate,length_estimate]=GMTPHD_estimation(weights_u, means_u,t_ini_u,length_u,Lscan,means_k_old_u)

%This estimator takes the mean of the trajectories corresponding to the highest peaks of the PHD.
%N_detections: number of detected trajectories

% X_estimate,t_b_estimate,length_estimate contain the trajectories
% estimates, their time of births and their lengths

%Author: Angel F. Garcia-Fernandez

N_detections=round(sum(weights_u));

X_estimate=cell(N_detections,1);
t_b_estimate=zeros(1,N_detections);
length_estimate=zeros(1,N_detections);
if(N_detections>0)   
    weights_u_copy=weights_u;
    for i=1:N_detections
        [~,index]=max(weights_u_copy);
        weights_u_copy(index)=0; %We set the chosen value to zero
        
        if(length_u(index)<=Lscan)
            X_estimate{i}=means_u(1:4*length_u(index),index);
        else
            X_estimate{i}=[ means_k_old_u(4*t_ini_u(index)-3:4*t_ini_u(index)-3+4*length_u(index)-4*Lscan-1,index);means_u(:,index)];
        end
    
        t_b_estimate(i)=t_ini_u(index);
        length_estimate(i)=length_u(index);        
    end    
end
