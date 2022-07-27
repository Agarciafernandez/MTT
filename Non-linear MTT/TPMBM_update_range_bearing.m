function filter_upd=TPMBM_update_range_bearing(filter_pred,z,p_d,k,gating_threshold,intensity_clutter,Nhyp_max,Lscan,kappa,Ap_kappa,x_s,weights_sp,R_kf,Nx,Nz,Nit_iplf,log_denom_likelihood,threshold_iplf,log_denom_Gauss,mesurement_dir_flag,likelihood_corr_flag)

%TPMBM update for alive trajectories with range-bearings measurements 
%IPLF implementation with optional likelihood improvement

%Author: Angel F. Garcia-Fernandez

Nprev_tracks=length(filter_pred.tracks);

%Poisson update (all trajectories are alive). Only the weight changes
Ncom_Pois=length(filter_pred.Pois);
filter_upd.Pois=filter_pred.Pois;

%change (with non-constant p_d)
for i=1:Ncom_Pois
    %We take the mean for the current target
    mean_i= filter_pred.Pois{i}.meanPois(end-Nx+1:end);
    %We calculate the average probability of detection with each method
    if(likelihood_corr_flag)
        %We take the covariance for the current target
        cov_i=filter_pred.Pois{i}.covPois(end-Nx+1:end,end-Nx+1:end);
        p_d_mean=integral_pd(mean_i,cov_i,p_d,weights_sp,Nx,x_s);
    else
        %We take the value at the mean
        mean_dist=sqrt((mean_i(1)-x_s(1)).^2+(mean_i(3)-x_s(2)).^2);
        p_d_mean=p_d(mean_dist);
    end
    filter_upd.Pois{i}.weightPois=(1-p_d_mean)*filter_pred.Pois{i}.weightPois;
end


%New track generation
Nnew_tracks=size(z,2); %Number of new tracks (some may have existence probability equal to zero)


for m=1:size(z,2)
    z_m=z(:,m);
    
    %We go through all Poisson components
    
    indices_new_tracks=[];
    for i=1:length(filter_pred.Pois)
        
        %Mean and covariance of the last target state
        mean_i=filter_pred.Pois{i}.meanPois(end-Nx+1:end);
        cov_i=filter_pred.Pois{i}.covPois(end-Nx+1:end,end-Nx+1:end);
        
        [z_pred_i,S_pred_i]=Predicted_range_bearings_VM(mean_i,cov_i,weights_sp,weights_sp(1),Nx,Nz,kappa,Ap_kappa,x_s);
        S_pred_i=S_pred_i+R_kf;
        
        
        maha=(z_m-z_pred_i)'/S_pred_i*(z_m-z_pred_i);
        
        if(maha<gating_threshold)
            %A track is created
            indices_new_tracks=[indices_new_tracks,i];
        end
        
    end
    
    if(~isempty(indices_new_tracks))
        %A new track is created for this measurement, based on the
        %component with highest weight (absorption)
        weightB=0;
        weight_i_max=0;
       
        
        for i=1:length(indices_new_tracks)
            index=indices_new_tracks(i);
            %Last target mean and covariance matrix
            mean_i=filter_pred.Pois{index}.meanPois(end-Nx+1:end);
            cov_i=filter_pred.Pois{index}.covPois(end-Nx+1:end,end-Nx+1:end);
            
            if(mesurement_dir_flag==1)
                %Measurement directed initialisation
                [mean_ini,P_ini]=Measurement_directed_ini_range_bearing(x_s,z_m,kappa,weights_sp(1));
            else
                mean_ini=mean_i;
                P_ini=cov_i;
            end
            
            [x_u_i,P_u_i,likelihood_i]=PL_range_bearing_update(mean_i, cov_i,mean_ini,P_ini,z_m,kappa,Ap_kappa,x_s, weights_sp,Nx,Nz,R_kf,Nit_iplf,...
                likelihood_corr_flag,p_d,log_denom_likelihood,threshold_iplf,log_denom_Gauss);
            
            
            weight_i=likelihood_i*filter_pred.Pois{index}.weightPois;
            weightB=weightB+weight_i;
            
            if(weight_i>=weight_i_max)
                index_i_max=index;
                x_u_i_max=x_u_i;
                P_u_i_max=P_u_i;
                weight_i_max=weight_i;
            end
            
        end
        
        %Update of the whole trajctory (or the last Lscan time steps)
        
        length_k_i=filter_pred.Pois{index_i_max}.length_Pois;
        min_length=min([Lscan,length_k_i]);
        
        
        if(min_length>1)
            %We use the partitioned update to update past states of the
            %trajectoreis
            %Prior mean and covariance
            mean_i=filter_pred.Pois{index_i_max}.meanPois(end-Nx*min_length+1:end);
            cov_i=filter_pred.Pois{index_i_max}.covPois(end-Nx*min_length+1:end,end-Nx*min_length+1:end);
            
            [mean_u,covB]=Update_past_trajectory(mean_i,cov_i,x_u_i_max,P_u_i_max,Nx);
            
        else
            mean_u=x_u_i_max;
            covB=P_u_i_max;
        end
        
        meanB=filter_pred.Pois{index_i_max}.meanPois;
        meanB(end-Nx*min_length+1:end)=mean_u;
        
        
        eB=weightB;
        
        %We add clutter to weight
        weightB=weightB+intensity_clutter;
        eB=eB/weightB;
        
        
        filter_upd.tracks{Nprev_tracks+m}.meanB{1}=meanB;
        filter_upd.tracks{Nprev_tracks+m}.covB{1}=covB;
        filter_upd.tracks{Nprev_tracks+m}.eB=eB;
        filter_upd.tracks{Nprev_tracks+m}.t_ini=k;
        filter_upd.tracks{Nprev_tracks+m}.aHis{1}=m;
        filter_upd.tracks{Nprev_tracks+m}.weightBLog=log(weightB);
        filter_upd.tracks{Nprev_tracks+m}.t_b=filter_pred.Pois{index_i_max}.t_bPois;
        filter_upd.tracks{Nprev_tracks+m}.length=filter_pred.Pois{index_i_max}.length_Pois;
        
        
    else
        %The Bernoulli component has existence probability zero (it is
        %clutter). It will be removed by pruning
        weightB=intensity_clutter;
        
        filter_upd.tracks{Nprev_tracks+m}.eB=0;
        filter_upd.tracks{Nprev_tracks+m}.weightBLog=log(weightB);
        
        filter_upd.tracks{Nprev_tracks+m}.meanB{1}=zeros(Nx,1);
        filter_upd.tracks{Nprev_tracks+m}.covB{1}=zeros(Nx,Nx);
        filter_upd.tracks{Nprev_tracks+m}.t_ini=k;
        filter_upd.tracks{Nprev_tracks+m}.aHis{1}=m;
        
    end
    
end




%We update all previous tracks with the measurements

for i=1:Nprev_tracks
    Nhyp_i=length(filter_pred.tracks{i}.eB);
    filter_upd.tracks{i}.t_ini=filter_pred.tracks{i}.t_ini;
    
    %We copy the birth time and length of the trajectory
    filter_upd.tracks{i}.t_b=filter_pred.tracks{i}.t_b;
    filter_upd.tracks{i}.length=filter_pred.tracks{i}.length;
    
    for j=1:Nhyp_i %We go through all hypotheses
        
        mean_j_traj=filter_pred.tracks{i}.meanB{j};
        cov_j_traj=filter_pred.tracks{i}.covB{j};
        eB_j=filter_pred.tracks{i}.eB(j);
        aHis_j=filter_pred.tracks{i}.aHis{j};
        
        
        %Misdetection hypotheses  Mean and covariance remain unchanged. Only the weights
        %change (with non-constant p_d)
        filter_upd.tracks{i}.meanB{j}=mean_j_traj;
        filter_upd.tracks{i}.covB{j}=cov_j_traj;
        
        %We take the mean and covariance of the last target state
        mean_j=mean_j_traj(end-Nx+1:end);
        cov_j=cov_j_traj(end-Nx+1:end,end-Nx+1:end);
        if(likelihood_corr_flag)
            p_d_mean=integral_pd(mean_j,cov_j,p_d,weights_sp,Nx,x_s);
        else
            %We take the value at the mean
            mean_dist=sqrt((mean_j(1)-x_s(1)).^2+(mean_j(3)-x_s(2)).^2);
            p_d_mean=p_d(mean_dist);
        end
        
        
        filter_upd.tracks{i}.eB(j)=eB_j*(1-p_d_mean)/(1-eB_j+eB_j*(1-p_d_mean));
        filter_upd.tracks{i}.aHis{j}=[aHis_j,0];
        
        filter_upd.tracks{i}.weightBLog_k(j)=log(1-eB_j+eB_j*(1-p_d_mean)); %Weight only at time step k (removing previous weight)
        
        
        %Predicted measurement and covariance matrix
        [z_pred_j,S_pred_j]=Predicted_range_bearings_VM(mean_j,cov_j,weights_sp,weights_sp(1),Nx,Nz,kappa,Ap_kappa,x_s);
        S_pred_j=S_pred_j+R_kf;      
        
        inv_S_pred_j=inv(S_pred_j);
                   
        length_k_i=filter_pred.tracks{i}.length;
        min_length=min([Lscan,length_k_i]);
   
        %We go through all measurements
        for m=1:size(z,2)
            index_hyp=j+Nhyp_i*m;
            z_m=z(:,m);
            
            
            maha=(z_m-z_pred_j)'*inv_S_pred_j*(z_m-z_pred_j);
            
            if(maha<gating_threshold)
                %We update the current target state with this measurement
                [x_u_c,P_u_c,likelihood_i]=PL_range_bearing_update(mean_j, cov_j,mean_j,cov_j,z_m,kappa,Ap_kappa,x_s, weights_sp,Nx,Nz,R_kf,Nit_iplf,likelihood_corr_flag,p_d,...
                    log_denom_likelihood,threshold_iplf,log_denom_Gauss);
                
                %We update past states of the trajectory (within the
                %L-scan window)
                if(min_length>1)                  
                    [x_u,P_u]=Update_past_trajectory(mean_j_traj(end-Nx*min_length+1:end),cov_j_traj,x_u_c,P_u_c,Nx);
                else
                    x_u=x_u_c;
                    P_u=P_u_c;
                end
     
                %We have to account for the values outslide the Lscan
                %window
                x_u_traj=mean_j_traj;
                x_u_traj(end-Nx*min_length+1:end)=x_u;
                
                filter_upd.tracks{i}.meanB{index_hyp}=x_u_traj;
                filter_upd.tracks{i}.covB{index_hyp}=P_u;
                filter_upd.tracks{i}.eB(index_hyp)=1;
                filter_upd.tracks{i}.aHis{index_hyp}=[aHis_j,m];
                
                %Update of the weights
                filter_upd.tracks{i}.weightBLog_k(index_hyp)=log(eB_j*likelihood_i) ;
                
                
            else
                filter_upd.tracks{i}.weightBLog_k(index_hyp)=-Inf;
                
                
            end
            
        end
        
    end
    
end

%Update of global hypotheses

if(Nprev_tracks==0)
    %No measurement hypothesised to be received from a target yet (No
    %previous tracks)
    %We create one global hypothesis
    filter_upd.globHypWeight=1;
    filter_upd.globHyp=ones(1,Nnew_tracks);
    
    if(Nnew_tracks==0)
        filter_upd.tracks=cell(0,1);
    end
    
else
    
    globWeightLog=[];
    globHyp=[];
    
    for p=1:length(filter_pred.globHypWeight)
        %We create cost matrix for each global hypothesis.
        cost_matrix_log=-Inf(Nprev_tracks+Nnew_tracks,size(z,2));
        cost_misdetection=zeros(Nprev_tracks,1);
        
        for i=1:Nprev_tracks
            index_hyp=filter_pred.globHyp(p,i); %Hypothesis for track i in p global hypothesis
            Nhyp_i=length(filter_pred.tracks{i}.eB);
            %We generate the cost matrix for measurements
            
            if(index_hyp~=0)
                
                index_max=length(filter_upd.tracks{i}.weightBLog_k);
                indices=index_hyp+Nhyp_i*(1:size(z,2));
                indices_c=indices<=index_max;
                indices_c=indices(indices_c);
                
                weights_log=filter_upd.tracks{i}.weightBLog_k(indices_c);
                
                %We remove the weight of the misdetection hypothesis to use Murty (Later this weight is added).
                cost_matrix_log(i,1:length(indices_c))=weights_log-filter_upd.tracks{i}.weightBLog_k(index_hyp);
                cost_misdetection(i)=filter_upd.tracks{i}.weightBLog_k(index_hyp);
                
                
            end
        end
        
        %New targets
        for i=Nprev_tracks+1:Nnew_tracks+Nprev_tracks
            weights_log=filter_upd.tracks{i}.weightBLog;
            index=filter_upd.tracks{i}.aHis{1};
            cost_matrix_log(i,index)=weights_log;
        end
        
        
        %We remove -Inf rows and columns for performing optimal assignment. We take them
        %into account for indexing later.
        %Columns that have only one value different from Inf are not fed
        %into Murty either.
        
        cos_matrix_inf=~isinf(cost_matrix_log);
        indices_stay_column=find(sum(cos_matrix_inf,1)>1); %There is one after inequality
        
        indices_stay_row=find(sum(cos_matrix_inf(:,indices_stay_column),2)>0);    %There is a zero after inequality
        
        fixed_indices_stay_column=find(sum(cos_matrix_inf,1)==1); %These belong to the final hypotheses but are not fed into Murty's algorithms
        [fixed_indices_stay_row,~]=find(cos_matrix_inf(:,fixed_indices_stay_column));
        
        
        %cost_fixed is the overall cost of the hypotheses that always belong to the output
        cost_fixed=sum(cost_matrix_log(fixed_indices_stay_row+size(cost_matrix_log,1)*(fixed_indices_stay_column'-1)));
        cost_matrix_log_trimmed=cost_matrix_log(indices_stay_row,indices_stay_column);
        
        
        
        globWeightLog_pred=log(filter_pred.globHypWeight(p));
        
        if(isempty(cost_matrix_log_trimmed))
            %One hypothesis: All targets are misdetected, no need for Murty
            opt_indices_trans=zeros(1,size(z,2));
            opt_indices_trans(fixed_indices_stay_column)=fixed_indices_stay_row;
            globWeightLog=[globWeightLog,sum(cost_misdetection)+cost_fixed+globWeightLog_pred];
            
        else
            %Number of new global hypotheses from this global hypothesis
            kbest=ceil(Nhyp_max*filter_pred.globHypWeight(p));
            
            
            %Option 1: We run Murty algorithm (making use of Hungarian algorithm to
            %solve the assginment problem). We call the function by using transpose and negative value
            %[opt_indices,nlcost]=murty(-cost_matrix_log_trimmed',kbest);
            
            %Option 2: We run MURTY algorithm (making use of Hungarian algorithm to
            %solve the assginment problem). We call the function by using transpose and negative value
            %kBest2DAssign by David F. Crouse
            [opt_indices,~,nlcost]=kBest2DAssign(-cost_matrix_log_trimmed',kbest);
            opt_indices=opt_indices';
            nlcost=nlcost';
            
            %Optimal indices without removing Inf rows
            opt_indices_trans=Inf(size(opt_indices,1),size(z,2));
            for i=1:size(opt_indices,1)
                opt_indices_trans(i,indices_stay_column)=indices_stay_row(opt_indices(i,:));
                
                %We add the single trajectory hypotheses that belong to all k
                %max
                opt_indices_trans(i,fixed_indices_stay_column)=fixed_indices_stay_row;
            end
            
            globWeightLog=[globWeightLog,-nlcost+sum(cost_misdetection)+cost_fixed+globWeightLog_pred];
        end
        
        %We write the corresponding global hypothesis
        globHypProv=zeros(size(opt_indices_trans,1),Nprev_tracks+Nnew_tracks);
        
        for i=1:Nprev_tracks
            index_hyp=filter_pred.globHyp(p,i); %Hypothesis for track i in p global hypothesis
            Nhyp_i=length(filter_pred.tracks{i}.eB);
            for j=1:size(opt_indices_trans,1)
                index_track= find(opt_indices_trans(j,:)==i);
                if(isempty(index_track))
                    index=index_hyp;
                else
                    index=index_hyp+Nhyp_i*index_track;
                end
                globHypProv(j,i)=index;
            end
        end
        
        
        for i=Nprev_tracks+1:Nnew_tracks+Nprev_tracks
            
            for j=1:size(opt_indices_trans,1)
                index_track=find(opt_indices_trans(j,:)==i);
                if(isempty(index_track))
                    index=0;
                else
                    index=1;
                end
                globHypProv(j,i)=index; %It is either 0 or 1
            end
        end
        
        globHyp=[globHyp;globHypProv];
        
        
    end
    filter_upd.globHyp=globHyp;
    %Normalisation of weights of global hypotheses
    globWeight=exp(globWeightLog-max(max(globWeightLog)));
    globWeight=globWeight/sum(globWeight);
    filter_upd.globHypWeight=globWeight;
    
end



