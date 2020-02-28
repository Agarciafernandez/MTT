function [weights_u, means_u,covs_u,Ncom_u,logical_active_u]=GMPHD_filter_update(weights_k, means_k,covs_k,Ncom_k,logical_active_k,z,H,R,p_d,l_clutter,Area)
%GMPHD filter update

%See "B.-N. Vo and W.-K. Ma, “The Gaussian mixture probability hypothesis
%density filter,” IEEE Transactions on Signal Processing, vol. 54, no. 11,pp. 4091–4104, Nov. 2006."


[Nz,Nx]=size(H);

index_active_k=find(logical_active_k==1);

z_pred_series=H*means_k(:,index_active_k);

S_pred_series=zeros(Nz,Nz,Ncom_k);
K_pred_series=zeros(Nx,Nz,Ncom_k);
P_u_series=zeros(Nx,Nx,Ncom_k);


for i=1:Ncom_k
    cov_i=covs_k(:,:,index_active_k(i));
    S_pred_i=H*cov_i*H'+R;
    S_pred_i=(S_pred_i+S_pred_i')/2;
    
    K_pred_i=cov_i*H'/S_pred_i;
    P_u_i=(eye(Nx)-K_pred_i*H)*cov_i;
    
    P_u_i=(P_u_i+P_u_i')/2;
    
    S_pred_series(:,:,i)=S_pred_i;
    K_pred_series(:,:,i)=K_pred_i;
    P_u_series(:,:,i)=P_u_i;
    
    
end

%Update (total number of terms is )
N_measurements=size(z,2);
Ncom_u=N_measurements*Ncom_k+Ncom_k;

weights_u=zeros(size(weights_k));
means_u=zeros(size(means_k));
covs_u=zeros(size(covs_k));
logical_active_u=false(size(logical_active_k));

logical_active_u(1:Ncom_u)=1;

%Update (non-detected components)
index_u=1;
for i=1:Ncom_k
    
   weight_i=(1-p_d)*weights_k(index_active_k(i));
   mean_i=means_k(:,index_active_k(i));
   cov_i=covs_k(:,:,index_active_k(i));
   
   weights_u(index_u)=weight_i;
   means_u(:,index_u)=mean_i;
   covs_u(:,:,index_u)=cov_i;
   index_u=index_u+1;    
end

%Update (detected components)

for p=1:N_measurements
    
    sum_weights_p=0;
    for i=1:Ncom_k
       
        weight_ip=p_d*weights_k(index_active_k(i))*mvnpdf(z(:,p),z_pred_series(:,i),S_pred_series(:,:,i));
        mean_ip=means_k(:,index_active_k(i))+K_pred_series(:,:,i)*(z(:,p)-z_pred_series(:,i));
        cov_ip=P_u_series(:,:,i);     
        weights_u(index_u)=weight_ip;
        means_u(:,index_u)=mean_ip;
        covs_u(:,:,index_u)=cov_ip;
        index_u=index_u+1;
        sum_weights_p=sum_weights_p+weight_ip;
     
    end
    
    %Division of the weights of the last Ncom_k components 
    intensity_clutter_p=l_clutter/(Area(1)*Area(2)); %Clutter es uniforme en el area
    weights_u(index_u-Ncom_k:index_u-1)=weights_u(index_u-Ncom_k:index_u-1)/(intensity_clutter_p+sum_weights_p);   
end