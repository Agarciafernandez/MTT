function [X_estimate,t_b_estimate,length_estimate]=GMTCPHD_estimation(weights_u, means_u,t_ini_u,length_u,Lscan,means_k_old_u,cardinality_u)

%Estimation for the GM-TCPHD filter
%We take the N highest peaks (according to the maximum of the cardinality) and report the mean of the trajectories at
%these peaks

%Author: Angel F. Garcia-Fernandez

[~,index_card]=max(cardinality_u);
N_detections=index_card-1;
X_estimate=cell(N_detections,1);
t_b_estimate=zeros(1,N_detections);
length_estimate=zeros(1,N_detections);

if(N_detections>0)  
    for i=1:N_detections
        [~,index_max]=max(weights_u);       
        if(length_u(index_max)<=Lscan)
            X_estimate{i}=means_u(1:4*length_u(index_max),index_max);
        else           
            X_estimate{i}=[ means_k_old_u(4*t_ini_u(index_max)-3:4*t_ini_u(index_max)-3+4*length_u(index_max)-4*Lscan-1,index_max);means_u(:,index_max)];      
        end       
        t_b_estimate(i)=t_ini_u(index_max);
        length_estimate(i)=length_u(index_max);      
        weights_u(index_max)=0;
    end
    
end







