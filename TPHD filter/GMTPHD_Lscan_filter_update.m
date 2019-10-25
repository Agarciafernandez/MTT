function [weights_u, means_u,covs_u,t_ini_u,length_u,logical_actives_u,means_k_old_u]=GMTPHD_Lscan_filter_update(weights_k, means_k,covs_k,t_ini_k,length_k,Ncom_k,...
    logical_actives_k,z,H,R,p_d,l_clutter,Area,Lscan,means_k_old)

%Update of the GMTPHD filter (with L-scan window) and only considering alive
%trajectories
%No gating is performed though it can be added
%Author: Angel Garcia-Fernandez

means_k_old_u=zeros(size(means_k_old));

[Nz,Nx]=size(H);

z_pred_series=zeros(Nz,1,Ncom_k);
S_pred_series=zeros(Nz,Nz,Ncom_k);
K_pred_series=zeros(Nx*Lscan,Nz,Ncom_k);
P_u_series=zeros(Nx*Lscan,Nx*Lscan,Ncom_k);


%We go through all PHD components
for i=1:Ncom_k
    
    %We consider the right indices over the L-scan window (considering
    %Nx=4)
    length_k_i=length_k(i);  
    if(length_k_i<Lscan)
        index_actual_i=Nx*length_k_i-(Nx-1):Nx*length_k_i;
        index_trajectory_i=1:Nx*length_k_i;
    else
        index_actual_i=Nx*Lscan-(Nx-1):Nx*Lscan;
        index_trajectory_i=1:Nx*Lscan;       
    end
    
    %Predicted measurement
    z_pred_i=H*means_k(index_actual_i,i);    
    
    %We construct the measurement matrix for the trajectory    
    min_length=min([Lscan,length_k_i]);
    H_trajectory_i=[zeros(Nz,Nx*(min_length-1)),H];
    
    
        
    %Predicted covariance matrix of the measurement    
    cov_i=covs_k(index_trajectory_i,index_trajectory_i,i);
    S_pred_i=H_trajectory_i*cov_i*H_trajectory_i'+R;
    S_pred_i=(S_pred_i+S_pred_i')/2;
    
    %Update of the covariance matrix
    K_pred_i=cov_i*H_trajectory_i'/S_pred_i;
    P_u_i=(eye(length(index_trajectory_i))-K_pred_i*H_trajectory_i)*cov_i;
    
    P_u_i=(P_u_i+P_u_i')/2;
    
    %We store the output of the Kalman filter update parameters 
    z_pred_series(:,1,i)=z_pred_i;
    S_pred_series(:,:,i)=S_pred_i;
    K_pred_series(index_trajectory_i,:,i)=K_pred_i;
    P_u_series(index_trajectory_i,index_trajectory_i,i)=P_u_i;
    
    
end

%Ncom_u: Total number of PHD components after update
N_measurements=size(z,2);
Ncom_u=N_measurements*Ncom_k+Ncom_k;

%Updated parameters of the PHD
weights_u=zeros(size(weights_k));
means_u=zeros(size(means_k));
covs_u=zeros(size(covs_k));
t_ini_u=zeros(size(t_ini_k));
length_u=zeros(size(length_k));
logical_actives_u=false(size(logical_actives_k));
logical_actives_u(1:Ncom_u)=1;

%Update (non detected components)
index_u=1; 
for i=1:Ncom_k
   weight_i=(1-p_d)*weights_k(i);
   mean_i=means_k(:,i);
   cov_i=covs_k(:,:,i);
   
   t_ini_k_i=t_ini_k(i);
   length_k_i=length_k(i);
   
   weights_u(index_u)=weight_i;
   means_u(:,index_u)=mean_i;
   covs_u(:,:,index_u)=cov_i;
   t_ini_u(index_u)=t_ini_k_i;
   length_u(index_u)=length_k_i;
   means_k_old_u(:,index_u)=means_k_old(:,i);
   
   index_u=index_u+1;    
end

%Update detected components (for loop over measurements)

for p=1:N_measurements    
    sum_weights_p=0;
    for i=1:Ncom_k
        
        length_k_i=length_k(i);       
        if(length_k_i<Lscan)
            index_trajectory_i=1:Nx*length_k_i;
        else
            index_trajectory_i=1:Nx*Lscan;
        end
    
        t_ini_k_i=t_ini_k(i);
        
        %Evaluation of the weights  
        gaussian_ev=evaluate2DGaussianpdf(z(:,p),z_pred_series(:,i),S_pred_series(:,:,i));
        weight_ip=p_d*weights_k(i)*gaussian_ev;
     
        
        mean_ip=means_k(index_trajectory_i,i)+K_pred_series(index_trajectory_i,:,i)*(z(:,p)-z_pred_series(:,i));
        cov_ip=P_u_series(:,:,i);     
        weights_u(index_u)=weight_ip;
        means_u(index_trajectory_i,index_u)=mean_ip;
        covs_u(:,:,index_u)=cov_ip;
        t_ini_u(index_u)=t_ini_k_i;
        length_u(index_u)=length_k_i;
        means_k_old_u(:,index_u)=means_k_old(:,i);
        
        index_u=index_u+1;
        sum_weights_p=sum_weights_p+weight_ip;
     
    end
    
    %Division of the last Ncom_k weights
    intensity_clutter_p=l_clutter/(Area(1)*Area(2)); %Clutter is uniformly distributed
    weights_u(index_u-Ncom_k:index_u-1)=weights_u(index_u-Ncom_k:index_u-1)/(intensity_clutter_p+sum_weights_p);
        
end




