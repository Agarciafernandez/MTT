function [weights_u, means_u,covs_u,t_ini_u,length_u,Ncom_u,logical_actives_u,cardinality_u,means_k_old_u]=GMTCPHD_Lscan_filter_update(weights_k, means_k,covs_k,t_ini_k,length_k,Ncom_k,logical_actives_k,z,H,R,p_d,l_clutter,Area,cardinality_pred,Lscan,means_k_old)

%Update step of the Gaussian mixture TCPHD filter
%No gating is performed though it can be added
%Author: Angel F. Garcia-Fernandez


%means outside the L-scan window
means_k_old_u=zeros(size(means_k_old));
Ncardinality_max=length(cardinality_pred)-1;


[Nz,Nx]=size(H);


%We compute the Kalman filter moments for the update for each PHD component
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

weights_u=zeros(size(weights_k));
means_u=zeros(size(means_k));
covs_u=zeros(size(covs_k));
t_ini_u=zeros(size(t_ini_k));
length_u=zeros(size(length_k));
logical_actives_u=false(size(logical_actives_k));

logical_actives_u(1:Ncom_u)=1;

%Update non-detected components
index_u=1;
for i=1:Ncom_k
    mean_i=means_k(:,i);
    cov_i=covs_k(:,:,i);
    length_k_i=length_k(i);
    t_ini_k_i=t_ini_k(i);
    means_u(:,index_u)=mean_i;
    covs_u(:,:,index_u)=cov_i;
    t_ini_u(index_u)=t_ini_k_i;
    length_u(index_u)=length_k_i;
    means_k_old_u(:,index_u)=means_k_old(:,i);
    index_u=index_u+1;
end

%Update detected componenents
weight_z=zeros(Ncom_k,N_measurements);
weights_pred=weights_k(1:Ncom_k)';
indeces_u_detect=zeros(Ncom_k,N_measurements);
for p=1:N_measurements
    for i=1:Ncom_k
        weight_z(i,p)=evaluate2DGaussianpdf(z(:,p),z_pred_series(:,i),S_pred_series(:,:,i));      
        indeces_u_detect(i,p)=index_u;     
        length_k_i=length_k(i);
        if(length_k_i<Lscan)
            index_trajectory_i=1:Nx*length_k_i;
        else
            index_trajectory_i=1:Nx*Lscan;
        end
        
        t_ini_k_i=t_ini_k(i);    
        mean_ip=means_k(index_trajectory_i,i)+K_pred_series(index_trajectory_i,:,i)*(z(:,p)-z_pred_series(:,i));
        cov_ip=P_u_series(:,:,i);
        means_u(index_trajectory_i,index_u)=mean_ip;
        covs_u(:,:,index_u)=cov_ip;
        t_ini_u(index_u)=t_ini_k_i;
        length_u(index_u)=length_k_i;
        means_k_old_u(:,index_u)=means_k_old(:,i);
        index_u=index_u+1;
        
    end
end

%Input of the elementary symmetric functions (Note we consider uniform
%clutter in an area of size Area(1)*Area(2)
Lambda_upd=zeros(N_measurements,1);
for i=1:N_measurements
    Lambda_upd(i)=p_d*weights_pred'*weight_z(:,i)*(Area(1)*Area(2));
end

%We calculate the elementary symmetric functions for all elements
%of Lambda_upd and with one element removed
esfs=Compute_esf(Lambda_upd);
esfs_minus_one=zeros(N_measurements,N_measurements);
for i=1:N_measurements
    Lambda_upd_i=[Lambda_upd(1:i-1);Lambda_upd(i+1:end)];
    esfs_minus_one(:,i)=Compute_esf(Lambda_upd_i);
end

%We calculate all upsilons
upsilon0=zeros(Ncardinality_max+1,1);
upsilon1=zeros(Ncardinality_max+1,1);
upsilon1_all_minus_one=zeros(Ncardinality_max+1,N_measurements);
for n=0:Ncardinality_max
    %Calculate upsilon0
    for j=0:min(N_measurements,n)
        upsilon0(n+1)=upsilon0(n+1)+exp(-l_clutter+(N_measurements-j)*log(l_clutter)+sum(log(n-j+1:n))+(n-j)*log(1-p_d)-j*log(sum(weights_pred)))*esfs(j+1);
    end
    %Calculate upsilon1
    for j=0:min(N_measurements,n-1)
        upsilon1(n+1)=upsilon1(n+1)+exp(-l_clutter+(N_measurements-j)*log(l_clutter)+sum(log(n-j:n))+(n-(j+1))*log(1-p_d)-(j+1)*log(sum(weights_pred)))*esfs(j+1);
    end
    %Calculate upsilon1_all_minus_one
    for i=1:N_measurements
        for j=0:min((N_measurements-1),n-1)
            upsilon1_all_minus_one(n+1,i)=upsilon1_all_minus_one(n+1,i)+exp(-l_clutter+((N_measurements-1)-j)*log(l_clutter)+sum(log(n-j:n))+(n-(j+1))*log(1-p_d)-(j+1)*...
                log(sum(weights_pred)))*esfs_minus_one(j+1,i);          
        end
    end
end

%Updated cardinality distribution
cardinality_u=upsilon0.*cardinality_pred;
cardinality_u=cardinality_u/sum(cardinality_u);
%Update of the PHD weights
%Missed detections
weights_u(1:Ncom_k)=(upsilon1'*cardinality_pred)/(upsilon0'*cardinality_pred)*(1-p_d)*weights_pred;
%Measurements
for i=1:N_measurements
    weights_u(indeces_u_detect(:,i))=(upsilon1_all_minus_one(:,i)'*cardinality_pred)/(upsilon0'*cardinality_pred)*p_d.*weight_z(:,i)*(Area(1)*Area(2)).*weights_pred;
end

