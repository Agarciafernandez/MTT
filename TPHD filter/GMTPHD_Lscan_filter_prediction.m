function   [weights_k, means_k,covs_k,t_ini_k,length_k,Ncom_k,logical_activas_k,means_k_old]=GMTPHD_Lscan_filter_prediction(weights_u, means_u,covs_u,t_ini_u,length_u,Ncom_u,logical_activas_u,F,Q,p_s,Lscan,means_k_old_in,k)

%Prediction of the GMTPHD filter (with L-scan window) and only considering alive
%trajectories
%Author: Angel Garcia-Fernandez

means_k_old=means_k_old_in;

Ncom_k=Ncom_u;
logical_activas_k=logical_activas_u;

index_actives_u=find(logical_activas_u==1);

weights_k=zeros(size(weights_u));
means_k=zeros(size(means_u));
covs_k=zeros(size(covs_u));
t_ini_k=zeros(size(t_ini_u));
t_ini_k(logical_activas_k)=t_ini_u(logical_activas_k);
%We increment the length of (alive) trajectories by one
length_k=zeros(size(length_u));
length_k(logical_activas_k)=length_u(logical_activas_k)+1;

Nx=size(F,1);
Nx_1=Nx-1;

%We go through all PHD components
for i=1:Ncom_u
    mean_i=means_u(:,index_actives_u(i));
    cov_i=covs_u(:,:,index_actives_u(i));
    weight_i=weights_u(index_actives_u(i));
    length_k_i=length_k(index_actives_u(i));
    
    if(length_k_i<Lscan+1) 
        %Current trajectory length is shorter than L-scan window
        index_actual_i=Nx*length_k_i-Nx_1:Nx*length_k_i;
        index_prev_i=index_actual_i-Nx;
        
        mean_i_pred=mean_i;
        mean_i_pred(index_actual_i)=F*mean_i_pred(index_prev_i);
        
        %Covariace
        cov_i_ant=cov_i(1:index_actual_i(1)-1,1:index_actual_i(1)-1);
        F_mode=zeros(Nx,size(cov_i_ant,1));
        F_mode(end-Nx_1:end,end-Nx_1:end)=F;
        cross=F_mode*cov_i_ant;
        cov_i_pred=cov_i;
        cov_i_pred(index_actual_i,index_actual_i)=F*cov_i_pred(index_prev_i,index_prev_i)*F'+Q;    
        cov_i_pred(1:index_actual_i(1)-1,index_actual_i)=cross';
        cov_i_pred(index_actual_i,1:index_actual_i(1)-1)=cross;
        
    else
       %We adjust the mean and covariance to the considered L-scan window
       mean_i_pred=zeros(Nx*Lscan,1);
       mean_i_pred(end-(Nx-1):end)=F*mean_i(end-(Nx-1):end);
       mean_i_pred(1:end-Nx)=mean_i(Nx+1:end);    
       means_k_old(Nx*k-Nx*Lscan+1:Nx*k-Nx*Lscan+Nx,index_actives_u(i))=mean_i(1:Nx);
       
       %Covariance   
       cov_i_ant=cov_i(Nx+1:end,Nx+1:end);     
       cov_i_pred=zeros(Nx*Lscan,Nx*Lscan);
       cov_i_pred(end-Nx_1:end,end-Nx_1:end)=F*cov_i(end-Nx_1:end,end-Nx_1:end)*F'+Q;
       
       if(Lscan>1)
           F_mode=zeros(Nx,Nx*(Lscan-1));
           F_mode(:,end-Nx_1:end)=F;
           cross=F_mode*cov_i_ant;
           cov_i_pred(end-Nx_1:end,1:end-Nx)=cross;
           cov_i_pred(1:end-Nx,end-Nx_1:end)=cross';
           cov_i_pred(1:end-Nx,1:end-Nx)=cov_i(Nx+1:end,Nx+1:end);
       end
             
    end
    

    weight_i_pred=weight_i*p_s;    
    weights_k(index_actives_u(i))=weight_i_pred;
    means_k(:,index_actives_u(i))=mean_i_pred;
    covs_k(:,:,index_actives_u(i))=cov_i_pred;
   
end
