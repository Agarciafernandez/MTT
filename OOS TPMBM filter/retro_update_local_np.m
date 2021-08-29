function filter_upd=retro_update_local_np(filter_retro,filter_upd,Nhyp_i,p_d,H,R,z,i,gating_threshold,Nx,Nz)

%Local hypotheses update for OOS measurement update. Case
%non-present/present, when the time of birth corresponds to ko

%Author: Angel Garcia-Fernandez

p_s2=filter_retro.p_s2;


for j=1:Nhyp_i %We go through all hypotheses
    
    mean_j=filter_retro.tracks{i}.meanB{j};
    cov_j=filter_retro.tracks{i}.covB{j};
    eB_j=filter_retro.tracks{i}.eB(j);
    aHis_j=filter_retro.tracks{i}.aHis{j};
    prob_length_j=filter_retro.tracks{i}.prob_length{j};
    %is_oos_j=filter_retro.tracks{i}.is_oos(:,j);
    
    
    mean_j2=filter_retro.tracks{i}.meanB2{j};
    cov_j2=filter_retro.tracks{i}.covB2{j};
    
    
    %We compute some required terms
    cross_p_d=p_d*p_s2;
    
    
    
    %Misdetection hypotheses
    filter_upd.tracks{i}.meanB{j}=mean_j;
    filter_upd.tracks{i}.covB{j}=cov_j;
    filter_upd.tracks{i}.eB(j)=eB_j*(1-cross_p_d)/(1-eB_j+eB_j*(1-cross_p_d));
    filter_upd.tracks{i}.aHis{j}=[aHis_j,0];
    filter_upd.tracks{i}.weightBLog_k(j)=log(1-eB_j+eB_j*(1-cross_p_d)); %Weight only at time step k (removing previous weight)
    
    filter_upd.tracks{i}.meanB2{j}=mean_j2;
    filter_upd.tracks{i}.covB2{j}=cov_j2;
    
    
    %Update the dead components
    filter_upd.tracks{i}.mean_past{j}=filter_retro.tracks{i}.mean_past{j};
    filter_upd.tracks{i}.cov_past{j}=filter_retro.tracks{i}.cov_past{j};
    
    if(length(prob_length_j)>1)
        filter_upd.tracks{i}.mean_past2{j}=filter_retro.tracks{i}.mean_past2{j};
        filter_upd.tracks{i}.cov_past2{j}=filter_retro.tracks{i}.cov_past2{j};
    end
    
   
    
    %Update length probabilities (misdetection)
    %sum_prob_length_j=(1-p_s2)*sum(prob_length_j)+(1-p_d)*p_s2*sum(prob_length_j);
    sum_prob_length_j=(1-p_s2)+(1-p_d)*p_s2;

    filter_upd.tracks{i}.prob_length{j}=(1-p_s2)*prob_length_j/sum_prob_length_j;
    filter_upd.tracks{i}.prob_length2{j}=(1-p_d)*p_s2*prob_length_j/sum_prob_length_j;
    
    
    %Detection hypotheses
    %We only need to update meanj2 and cov_j2 as these have the state at
    %OOS time
    %KF moments
    S_pred_j=H*cov_j2(end-Nx+1:end,end-Nx+1:end)*H'+R;
    z_pred_j=H*mean_j2(end-Nx+1:end);
    
  
    
    %Inverse S_pred_j
    Vs= chol(S_pred_j);
    log_det_S_pred_j= 2*log(prod(diag(Vs)));
    inv_sqrt_S= inv(Vs);
    inv_S_pred_j= inv_sqrt_S*inv_sqrt_S';
    
    
    
    min_length=size(cov_j2,1)/Nx;
    H_trajectory=[zeros(Nz,Nx*(min_length-1)),H];
    K_pred_j=cov_j2*H_trajectory'*inv_S_pred_j;
    
    
    %We go through all measurements
    for m=1:size(z,2)
        
        
        
        index_hyp=j+Nhyp_i*m;
        z_m=z(:,m);
        
        maha=(z_m-z_pred_j)'*inv_S_pred_j*(z_m-z_pred_j);
        
        if(maha<gating_threshold)
            %We create the component
            
            P_u=(eye(Nx*min_length)-K_pred_j*H_trajectory)*cov_j2;
            P_u=(P_u+P_u')/2;
            x_u=mean_j2;
            x_u(end-Nx*min_length+1:end)=mean_j2(end-Nx*min_length+1:end)+K_pred_j*(z_m-z_pred_j);
            
            
            filter_upd.tracks{i}.meanB{index_hyp}=mean_j; %meanB and covB are not updated
            filter_upd.tracks{i}.covB{index_hyp}=cov_j;
            filter_upd.tracks{i}.meanB2{index_hyp}=x_u;
            filter_upd.tracks{i}.covB2{index_hyp}=P_u;
            filter_upd.tracks{i}.eB(index_hyp)=1;
            filter_upd.tracks{i}.aHis{index_hyp}=[aHis_j,m];
            
            %Update of the weights
            quad=-0.5*maha;
            p_z_log=quad-1/2*log_det_S_pred_j-Nz*log(2*pi)/2;
            
            % weight_partial_log=log(eB_j*cross_p_d)+p_z_log;
            
            weight_partial=zeros(size(prob_length_j));
            
            weight_partial(1)=p_d*exp(p_z_log)*p_s2*prob_length_j(1);
            
            Npast_means=length(filter_retro.tracks{i}.mean_past{j});
            
            
            for p=1:Npast_means
                
                %mean_past{j}{p} and cov_past{j}{p} do not have a state at OOS time
                
                mean_past2_jp=filter_retro.tracks{i}.mean_past2{j}{p};
                cov_past2_jp=filter_retro.tracks{i}.cov_past2{j}{p};
                
                [mean_past2_jp_u,cov_past2_jp_u,p_z_log_jp_u]=update_past_trajectory(mean_past2_jp,cov_past2_jp,H,R,z_m,Nz,Nx);
                
                
                filter_upd.tracks{i}.mean_past2{index_hyp}{p}=mean_past2_jp_u;
                filter_upd.tracks{i}.cov_past2{index_hyp}{p}=cov_past2_jp_u;
                
                weight_partial(p+1)=p_d*exp(p_z_log_jp_u)*p_s2*prob_length_j(p+1);
                
            end
            
            filter_upd.tracks{i}.weightBLog_k(index_hyp)=log(eB_j*sum(weight_partial));
            
            
            %Update the length probabilities and mean_past
            filter_upd.tracks{i}.prob_length{index_hyp}=0;
            
            filter_upd.tracks{i}.prob_length2{index_hyp}=weight_partial/sum(weight_partial);
            
            filter_upd.tracks{i}.mean_past{index_hyp}=cell(0,1);
            filter_upd.tracks{i}.cov_past{index_hyp}=cell(0,1);
            
        else
            filter_upd.tracks{i}.weightBLog_k(index_hyp)=-Inf;
            
            
        end
        
    end
    
    
end

end

function  [x_u,P_u,p_z_log]=update_past_trajectory(mean_j,cov_j,H,R,z_m,Nz,Nx)


%KF moments
S_pred_j=H*cov_j(end-Nx+1:end,end-Nx+1:end)*H'+R;
z_pred_j=H*mean_j(end-Nx+1:end);

%Inverse S_pred_j
Vs= chol(S_pred_j);
log_det_S_pred_j= 2*log(prod(diag(Vs)));
inv_sqrt_S= inv(Vs);
inv_S_pred_j= inv_sqrt_S*inv_sqrt_S';



min_length=size(cov_j,1)/Nx;
H_trajectory=[zeros(Nz,Nx*(min_length-1)),H];
K_pred_j=cov_j*H_trajectory'*inv_S_pred_j;

%Mahalanobis distance
maha=(z_m-z_pred_j)'*inv_S_pred_j*(z_m-z_pred_j);


P_u=(eye(Nx*min_length)-K_pred_j*H_trajectory)*cov_j;
P_u=(P_u+P_u')/2;
x_u=mean_j;
x_u(end-Nx*min_length+1:end)=mean_j(end-Nx*min_length+1:end)+K_pred_j*(z_m-z_pred_j);

%Update of the weights
quad=-0.5*maha;
p_z_log=quad-1/2*log_det_S_pred_j-Nz*log(2*pi)/2;


end


