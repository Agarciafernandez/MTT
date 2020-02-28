function [weights_u, means_u,covs_u,Ncom_u,logical_activas_u,cardinality_u]=GMCPHD_filter_update(weights_k, means_k,covs_k,Ncom_k,logical_activas_k,z,H,R,p_d,l_clutter,Area,cardinality_pred)


%Update of the GM-CPHD filter

Ncardinality_max=length(cardinality_pred)-1;

[Nz,Nx]=size(H);

index_activas_k=find(logical_activas_k==1);

z_pred_series=H*means_k(:,index_activas_k);

S_pred_series=zeros(Nz,Nz,Ncom_k);
K_pred_series=zeros(Nx,Nz,Ncom_k);
P_u_series=zeros(Nx,Nx,Ncom_k);


for i=1:Ncom_k
    cov_i=covs_k(:,:,index_activas_k(i));
    S_pred_i=H*cov_i*H'+R;
    S_pred_i=(S_pred_i+S_pred_i')/2;    
    K_pred_i=cov_i*H'/S_pred_i;
    P_u_i=(eye(Nx)-K_pred_i*H)*cov_i;    
    P_u_i=(P_u_i+P_u_i')/2;  
    S_pred_series(:,:,i)=S_pred_i;
    K_pred_series(:,:,i)=K_pred_i;
    P_u_series(:,:,i)=P_u_i;
       
end



%Update (total number of terms is)
N_measurements=size(z,2);
Ncom_u=N_measurements*Ncom_k+Ncom_k;

weights_u=zeros(size(weights_k));
means_u=zeros(size(means_k));
covs_u=zeros(size(covs_k));
logical_activas_u=false(size(logical_activas_k));

logical_activas_u(1:Ncom_u)=1;

%Update (non detected PHD components)
index_u=1;
for i=1:Ncom_k
    mean_i=means_k(:,index_activas_k(i));
    cov_i=covs_k(:,:,index_activas_k(i));
    means_u(:,index_u)=mean_i;
    covs_u(:,:,index_u)=cov_i;
    index_u=index_u+1;
end

%Update (detected PHD components)
weight_z=zeros(Ncom_k,N_measurements);
weights_pred=weights_k(index_activas_k)';
indeces_u_detect=zeros(Ncom_k,N_measurements);
for p=1:N_measurements
    
    for i=1:Ncom_k
        
        weight_z(i,p)=mvnpdf(z(:,p),z_pred_series(:,i),S_pred_series(:,:,i));
        indeces_u_detect(i,p)=index_u;
        
        mean_ip=means_k(:,index_activas_k(i))+K_pred_series(:,:,i)*(z(:,p)-z_pred_series(:,i));
        cov_ip=P_u_series(:,:,i);
        means_u(:,index_u)=mean_ip;
        covs_u(:,:,index_u)=cov_ip;
        index_u=index_u+1;
        
    end
    
    
    
end

%Input of the elementary symmetric functions
Lambda_upd=zeros(N_measurements,1);
for i=1:N_measurements
    Lambda_upd(i)=p_d*weights_pred'*weight_z(:,i)*(Area(1)*Area(2));
end

%We calculate the elementary symmetric functions for all elements
%of Lambda_upd and with one element removed
esf_all=Compute_esf(Lambda_upd);
esf_all_minus_one=zeros(N_measurements,N_measurements);
for i=1:N_measurements
    Lambda_upd_i=[Lambda_upd(1:i-1);Lambda_upd(i+1:end)];
    esf_all_minus_one(:,i)=Compute_esf(Lambda_upd_i);
end

%We calculate all upsilons
upsilon0=zeros(Ncardinality_max+1,1);
upsilon1=zeros(Ncardinality_max+1,1);
upsilon1_all_minus_one=zeros(Ncardinality_max+1,N_measurements);
for n=0:Ncardinality_max
    %Calculate upsilon0
    for j=0:min(N_measurements,n)
        upsilon0(n+1)=upsilon0(n+1)+exp(-l_clutter+(N_measurements-j)*log(l_clutter)+sum(log(n-j+1:n))+(n-j)*log(1-p_d)-j*log(sum(weights_pred)))*esf_all(j+1);
    end
    %Calculate upsilon1
    for j=0:min(N_measurements,n-1)
        upsilon1(n+1)=upsilon1(n+1)+exp(-l_clutter+(N_measurements-j)*log(l_clutter)+sum(log(n-j:n))+(n-(j+1))*log(1-p_d)-(j+1)*log(sum(weights_pred)))*esf_all(j+1);
    end
    %Calculate upsilon1_all_minus_one
    for i=1:N_measurements
        for j=0:min((N_measurements-1),n-1)
            upsilon1_all_minus_one(n+1,i)=upsilon1_all_minus_one(n+1,i)+exp(-l_clutter+((N_measurements-1)-j)*log(l_clutter)+sum(log(n-j:n))+(n-(j+1))*log(1-p_d)-(j+1)*...
                log(sum(weights_pred)))*esf_all_minus_one(j+1,i);          
        end
    end
end

%Cardinality update
cardinality_u=upsilon0.*cardinality_pred;
cardinality_u=cardinality_u/sum(cardinality_u);
%Update of the PHD weights
%Missed detections
weights_u(1:Ncom_k)=(upsilon1'*cardinality_pred)/(upsilon0'*cardinality_pred)*(1-p_d)*weights_pred;
%Measurements
for i=1:N_measurements
    weights_u(indeces_u_detect(:,i))=(upsilon1_all_minus_one(:,i)'*cardinality_pred)/(upsilon0'*cardinality_pred)*p_d.*weight_z(:,i)*(Area(1)*Area(2)).*weights_pred;
end

