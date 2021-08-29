function    filter_upd=retro_update_local_before_ko(filter_retro,filter_upd,Nhyp_i,p_d,H,R,z,i,gating_threshold,Nx,Nz,ko,Lscan)


%Local hypotheses update for OOS measurement update. Case
%non-present/present, when the time of birth corresponds to ko

%Author: Angel Garcia-Fernandez

p_s1=filter_retro.p_s1;

t_bi=filter_upd.tracks{i}.t_b;

%We store t_f (t_final for the main trajectory) to make it easier to perform
%the marginalisation (done after the update)

filter_upd.tracks{i}.t_f=zeros(Nhyp_i*(size(z,2)+1),1);

for j=1:Nhyp_i %We go through all hypotheses
    mean_j=filter_retro.tracks{i}.meanB{j};
    cov_j=filter_retro.tracks{i}.covB{j};
    eB_j=filter_retro.tracks{i}.eB(j);
    aHis_j=filter_retro.tracks{i}.aHis{j};
    prob_length_j=filter_retro.tracks{i}.prob_length{j};
    
        
    t_f=filter_retro.tracks{i}.t_f(j);
    
    
    
    filter_upd.tracks{i}.t_f(j)=t_f;
    
    
    if(t_f<ko-1)
        %Trajectory is not alive at OOS time
        
        %Misdetection hypotheses
        filter_upd.tracks{i}.meanB{j}=mean_j;
        filter_upd.tracks{i}.covB{j}=cov_j;
        filter_upd.tracks{i}.eB(j)=eB_j;
        filter_upd.tracks{i}.aHis{j}=[aHis_j,0];
        filter_upd.tracks{i}.weightBLog_k(j)=0; %Weight only at time step k (removing previous weight)
        filter_upd.tracks{i}.prob_length{j}=prob_length_j;
        
        %Update of the dead components
        filter_upd.tracks{i}.mean_past{j}=filter_retro.tracks{i}.mean_past{j};
        filter_upd.tracks{i}.cov_past{j}=filter_retro.tracks{i}.cov_past{j};
        
        
        %Detection hypotheses
        %We go through all measurements
        for m=1:size(z,2)
            index_hyp=j+Nhyp_i*m;
            filter_upd.tracks{i}.weightBLog_k(index_hyp)=-Inf;
        end
    elseif(t_f==ko-1)
        %We have a case of present/non-present
        
        mean_j2=filter_retro.tracks{i}.meanB2{j};
        cov_j2=filter_retro.tracks{i}.covB2{j};
        
        %Effective p_d
        cross_p_d=p_d*p_s1*prob_length_j(1);
        
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
                

        %Update length probabilities misdetection
        prob_length_j_upd=prob_length_j;
        prob_length_j_upd(1)=(1-p_s1)*prob_length_j(1);
        prob_length_j2_upd=(1-p_d)*p_s1*prob_length_j(1);
        
        norm=sum(prob_length_j_upd)+prob_length_j2_upd;
        prob_length_j_upd=prob_length_j_upd/norm;
        prob_length_j2_upd=prob_length_j2_upd/norm;
        
        
        %Update length probabilities (misdetection)
  
        
        filter_upd.tracks{i}.prob_length{j}=prob_length_j_upd;
        filter_upd.tracks{i}.prob_length2{j}=prob_length_j2_upd;
        
        
        
        %Detection hypothesis
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
                
                filter_upd.tracks{i}.t_f(index_hyp)=t_f;
                
                %Update of the weights
                quad=-0.5*maha;
                p_z_log=quad-1/2*log_det_S_pred_j-Nz*log(2*pi)/2;
                
                filter_upd.tracks{i}.weightBLog_k(index_hyp)=log(eB_j*cross_p_d)+p_z_log;
                
                %Update the length probabilities and mean_past
                filter_upd.tracks{i}.prob_length{index_hyp}=0;
                
                filter_upd.tracks{i}.prob_length2{index_hyp}=1;
                
                filter_upd.tracks{i}.mean_past{index_hyp}=cell(0,1);
                filter_upd.tracks{i}.cov_past{index_hyp}=cell(0,1);
            else
                filter_upd.tracks{i}.weightBLog_k(index_hyp)=-Inf;
            end
        end
        
    elseif(t_f>=ko)
        %We have a case of present-present (though past
        %trajectories may be present/non-present)
        
        
        %In present/present, we do not obtain a mixture in the
        %retrodiction. We do the update on mean_j, cov_j
        
        
        %In the next for loop, we sum over all alphas (see paper) with the aim of
        %obtaining the effective probability of detection
        sum_alpha=prob_length_j(1);
        
        prob_length_j_upd=zeros(size(prob_length_j));
        prob_length2_j_upd=zeros(size(prob_length_j));
        
        prob_length_j_upd(1)=(1-p_d)*prob_length_j(1);
        
        %We continue with past trajectories
        Npast_means=length(filter_retro.tracks{i}.mean_past{j});
        
        for p=1:Npast_means
            
            t_f_p=filter_retro.tracks{i}.t_f_past{j}{p};
            
            
            filter_upd.tracks{i}.t_f_past{j}{p}=t_f_p; %We store this variable to make marginalisation easier.
            
            
            if(t_f_p==ko-1)
                %We have a case of present/non-present
                sum_alpha=sum_alpha+p_s1*prob_length_j(p+1);
                prob_length2_j_upd(p+1)=(1-p_d)*p_s1*prob_length_j(p+1);
                prob_length_j_upd(p+1)=(1-p_s1)*prob_length_j(p+1);
                
            elseif(t_f_p>=ko)
                %We have a case of present/present
                sum_alpha=sum_alpha+prob_length_j(p+1);
                prob_length_j_upd(p+1)=(1-p_d)*prob_length_j(p+1);
                
            else
                prob_length_j_upd(p+1)=prob_length_j(p+1);
                
            end
            
        end
        cross_p_d=p_d*sum_alpha;
        
        %To normalise prob_length and prob_length2
        norm_prob_length=sum(prob_length_j_upd)+sum(prob_length2_j_upd);
        
        
        %Misdetection hypothesis
        filter_upd.tracks{i}.meanB{j}=mean_j;
        filter_upd.tracks{i}.covB{j}=cov_j;
        filter_upd.tracks{i}.eB(j)=eB_j*(1-cross_p_d)/(1-eB_j+eB_j*(1-cross_p_d));
        filter_upd.tracks{i}.aHis{j}=[aHis_j,0];
        filter_upd.tracks{i}.weightBLog_k(j)=log(1-eB_j+eB_j*(1-cross_p_d)); %Weight only at time step k (removing previous weight)
        
        filter_upd.tracks{i}.mis{j}=1; %We flag this is a misdetection hypothesis. To make marginalisation easier.
        
        %Update the dead components
        filter_upd.tracks{i}.mean_past{j}=filter_retro.tracks{i}.mean_past{j};
        filter_upd.tracks{i}.cov_past{j}=filter_retro.tracks{i}.cov_past{j};
        
        if(isfield(filter_retro.tracks{i},'mean_past2') && length(filter_retro.tracks{i}.mean_past2)>=j)
            filter_upd.tracks{i}.mean_past2{j}=filter_retro.tracks{i}.mean_past2{j};
            filter_upd.tracks{i}.cov_past2{j}=filter_retro.tracks{i}.cov_past2{j};
            prob_length2_j_upd=prob_length2_j_upd/norm_prob_length;
            filter_upd.tracks{i}.prob_length2{j}=prob_length2_j_upd;
        end
        
        
        %Update length probabilities (misdetection)
        prob_length_j_upd=prob_length_j_upd/norm_prob_length;
        
        filter_upd.tracks{i}.prob_length{j}=prob_length_j_upd;
        
        %Detection hypotheses
        %For the main trajectory, we only need to update meanj and cov_j as these have the state at
        %OOS time
        
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
        %We go through all measurements
        for m=1:size(z,2)
            
            
            index_hyp=j+Nhyp_i*m;
            z_m=z(:,m);
            
            maha=(z_m-z_pred_j)'*inv_S_pred_j*(z_m-z_pred_j);
            
            if(maha<gating_threshold)
                %We create the component
                
                P_u=(eye(Nx*min_length)-K_pred_j*H_trajectory)*cov_j;
                P_u=(P_u+P_u')/2;
                x_u=mean_j;
                x_u(end-Nx*min_length+1:end)=mean_j(end-Nx*min_length+1:end)+K_pred_j*(z_m-z_pred_j);
                
                filter_upd.tracks{i}.meanB{index_hyp}=x_u; %Update of x_u and P_u
                filter_upd.tracks{i}.covB{index_hyp}=P_u;
                
                filter_upd.tracks{i}.eB(index_hyp)=1;
                filter_upd.tracks{i}.aHis{index_hyp}=[aHis_j,m];
                
                filter_upd.tracks{i}.t_f(index_hyp)=t_f;
                
                filter_upd.tracks{i}.mis{index_hyp}=0;
                
                %Update of the weights
                quad=-0.5*maha;
                p_z_log=quad-1/2*log_det_S_pred_j-Nz*log(2*pi)/2;
                
                weight_partial=zeros(size(prob_length_j));
                weight_partial2=zeros(size(prob_length_j)); %This is for the alternative hypotheses (when one trajectory can generate two trajectories
                
                weight_partial(1)=p_d*exp(p_z_log)*prob_length_j(1);
                
                if(Npast_means==0 || Lscan<=2)
                    filter_upd.tracks{i}.mean_past{index_hyp}=cell(0,1);
                    filter_upd.tracks{i}.cov_past{index_hyp}=cell(0,1);
                    
                    
                else
                    
                    
                    for p=1:min(Npast_means,Lscan-2)
                        t_f_p=t_bi+length(filter_retro.tracks{i}.mean_past{j}{p})/Nx-1;
                        
                        filter_upd.tracks{i}.t_f_past{index_hyp}{p}=t_f_p;
                        
                        
                        if(t_f_p==ko-1)
                            %We have a case of present/non-present
                            %The one without augmentation passes unchanged, but
                            %with weight zero
                            
                            %Trajectory without state augmentation
                            %weight_partial(p+1)=0; No need to anything as this
                            %weight is already set to zero
                            mean_past_jp=filter_retro.tracks{i}.mean_past{j}{p};
                            cov_past_jp=filter_retro.tracks{i}.cov_past{j}{p};
                            filter_upd.tracks{i}.mean_past{index_hyp}{p}=mean_past_jp;
                            filter_upd.tracks{i}.cov_past{index_hyp}{p}=cov_past_jp;
                            
                            %Trajectory with state augmentation
                            mean_past2_jp=filter_retro.tracks{i}.mean_past2{j}{p};
                            cov_past2_jp=filter_retro.tracks{i}.cov_past2{j}{p};
                            
                            [mean_past2_jp_u,cov_past2_jp_u,p_z_log_jp_u]=update_past_trajectory(mean_past2_jp,cov_past2_jp,H,R,z_m,Nz,Nx);
                            
                            
                            filter_upd.tracks{i}.mean_past2{index_hyp}{p}=mean_past2_jp_u;
                            filter_upd.tracks{i}.cov_past2{index_hyp}{p}=cov_past2_jp_u;
                            
                            weight_partial2(p+1)=p_d*exp(p_z_log_jp_u)*p_s1*prob_length_j(p+1);
                            
                        elseif(t_f_p>=ko)
                            %Trajectory is present at both time steps
                            
                            mean_past_jp=filter_retro.tracks{i}.mean_past{j}{p};
                            cov_past_jp=filter_retro.tracks{i}.cov_past{j}{p};
                            
                            [mean_past_jp_u,cov_past_jp_u,p_z_log_jp_u]=update_past_trajectory(mean_past_jp,cov_past_jp,H,R,z_m,Nz,Nx);
                            
                            
                            filter_upd.tracks{i}.mean_past{index_hyp}{p}=mean_past_jp_u;
                            filter_upd.tracks{i}.cov_past{index_hyp}{p}=cov_past_jp_u;
                            
                            weight_partial(p+1)=p_d*exp(p_z_log_jp_u)*prob_length_j(p+1);
                        else
                            %Trajectory {j}{p} died before ko-1
                            %The trajectory does not change
                            
                            mean_past_jp=filter_retro.tracks{i}.mean_past{j}{p};
                            cov_past_jp=filter_retro.tracks{i}.cov_past{j}{p};
                            filter_upd.tracks{i}.mean_past{index_hyp}{p}=mean_past_jp;
                            filter_upd.tracks{i}.cov_past{index_hyp}{p}=cov_past_jp;
                            
                            weight_partial(p+1)=0; %as this is a detection hypothesis
                            
                        end
                        
                    end
                end
                
                filter_upd.tracks{i}.weightBLog_k(index_hyp)=log(eB_j*sum(weight_partial+weight_partial2));
                
                %Update the length probabilities and mean_past
                filter_upd.tracks{i}.prob_length{index_hyp}=weight_partial/sum(weight_partial+weight_partial2);
                
                filter_upd.tracks{i}.prob_length2{index_hyp}=weight_partial2/sum(weight_partial+weight_partial2);
                
                
            else
                filter_upd.tracks{i}.weightBLog_k(index_hyp)=-Inf;
                filter_upd.tracks{i}.mean_past{index_hyp}=cell(0,1);
                filter_upd.tracks{i}.cov_past{index_hyp}=cell(0,1);
                
            end
            
            
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

maha=(z_m-z_pred_j)'*inv_S_pred_j*(z_m-z_pred_j);


P_u=(eye(Nx*min_length)-K_pred_j*H_trajectory)*cov_j;
P_u=(P_u+P_u')/2;
x_u=mean_j;
x_u(end-Nx*min_length+1:end)=mean_j(end-Nx*min_length+1:end)+K_pred_j*(z_m-z_pred_j);

%Update of the weights
quad=-0.5*maha;
p_z_log=quad-1/2*log_det_S_pred_j-Nz*log(2*pi)/2;


end


