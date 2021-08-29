function filter_upd=TPMBM_all_retro_update(filter_retro,z,H,R,p_d,gating_threshold,intensity_clutter,Nhyp_max,Lscan)
%This function performs the TPMBM update with an OOS measurement, based on
%the retrodicted TPMBM density

%Author: Angel Garcia-Fernandez

ko=filter_retro.ko;

filter_upd.ko=ko;

[Nz,Nx]=size(H);
Nprev_tracks=length(filter_retro.tracks);

%Poisson update
Ncom_Pois=length(filter_retro.Pois);

filter_upd.Pois=filter_retro.Pois;

%We multiply by 1-p_d the weight of the trajectories that are alive at OOS time
for i=1:Ncom_Pois
    if(filter_retro.Pois{i}.is_oos)
        filter_upd.Pois{i}.weightPois=(1-p_d)*filter_upd.Pois{i}.weightPois;
    end
end




%New track generation
Nnew_tracks=size(z,2); %Number of new tracks (some may have existence probability equal to zero)


for m=1:size(z,2)
    z_m=z(:,m);
    
    %We go through all Poisson components
    
    indices_new_tracks=[];
    for i=1:length(filter_retro.Pois)
        
        %if(filter_retro.Pois{i}.is_oos &&
        %filter_retro.Pois{i}.weightPois>0) (This was used for testing)
        if(filter_retro.Pois{i}.is_oos)
            %We only consider PPP components with OOS state
            
            mean_i=filter_retro.Pois{i}.meanPois(end-Nx+1:end);
            cov_i=filter_retro.Pois{i}.covPois(end-Nx+1:end,end-Nx+1:end);
            
            
            S_pred_i=H*cov_i*H'+R;
            z_pred_i=H*mean_i;
            
            maha=(z_m-z_pred_i)'/S_pred_i*(z_m-z_pred_i);
            
            if(maha<gating_threshold)
                %A track is created
                indices_new_tracks=[indices_new_tracks,i];
            end
        end
        
    end
    
    if(~isempty(indices_new_tracks))
        %A new track is created for this measurement, based on the
        %component with highest weight (absorption)
        weightB=0;
        weight_i_max=0;
        z_pred_i_max=0;
        S_pred_i_max=zeros(Nz);
        
        for i=1:length(indices_new_tracks)
            %The trajectories in this for loop already have OOS measurement
            index=indices_new_tracks(i);
            mean_i=filter_retro.Pois{index}.meanPois(end-Nx+1:end);
            cov_i=filter_retro.Pois{index}.covPois(end-Nx+1:end,end-Nx+1:end);
            
            
            S_pred_i=H*cov_i*H'+R;
            z_pred_i=H*mean_i;
            
            weight_i=p_d*mvnpdf(z_m,z_pred_i,S_pred_i)*filter_retro.Pois{index}.weightPois;
            weightB=weightB+weight_i;
            
            if(weight_i>weight_i_max)
                index_i_max=index;
                z_pred_i_max=z_pred_i;
                S_pred_i_max=S_pred_i;
                weight_i_max=weight_i;
                
            end
            
        end
        
        %Update of the whole trajectory (or the last Lscan time steps)
        z_pred_i=z_pred_i_max;
        S_pred_i=S_pred_i_max;
        
        cov_i=filter_retro.Pois{index_i_max}.covPois;
        min_length=size(filter_retro.Pois{index_i_max}.covPois,1)/Nx;
        
        H_trajectory=[zeros(Nz,Nx*(min_length-1)),H];
        
        
        mean_i=filter_retro.Pois{index_i_max}.meanPois(end-Nx*min_length+1:end);
        
        K_pred_i=cov_i*H_trajectory'/S_pred_i;
        P_u_i=(eye(length(mean_i))-K_pred_i*H_trajectory)*cov_i;
        P_u_i=(P_u_i+P_u_i')/2;
        
        mean_u=mean_i+K_pred_i*(z_m-z_pred_i);
        
        
        meanB=filter_retro.Pois{index_i_max}.meanPois;
        meanB(end-Nx*min_length+1:end)=mean_u;
        covB=P_u_i;
        
        
        eB=weightB;
        
        %We add clutter to weight
        weightB=weightB+intensity_clutter;
        eB=eB/weightB;
        
        
        filter_upd.tracks{Nprev_tracks+m}.meanB{1}=meanB;
        filter_upd.tracks{Nprev_tracks+m}.covB{1}=covB;
        filter_upd.tracks{Nprev_tracks+m}.eB=eB;
        filter_upd.tracks{Nprev_tracks+m}.t_ini=-1; %To flag this is an OOS state
        filter_upd.tracks{Nprev_tracks+m}.aHis{1}=m;
        filter_upd.tracks{Nprev_tracks+m}.weightBLog=log(weightB);
        filter_upd.tracks{Nprev_tracks+m}.t_b=filter_retro.Pois{index_i_max}.t_bPois;
        filter_upd.tracks{Nprev_tracks+m}.length=filter_retro.Pois{index_i_max}.length_Pois;
        
        %When we track all trajectories, we also consider variable
        %prob_length
        filter_upd.tracks{Nprev_tracks+m}.prob_length{1}=1;
        filter_upd.tracks{Nprev_tracks+m}.mean_past{1}=cell(0,1);
        filter_upd.tracks{Nprev_tracks+m}.cov_past{1}=cell(0,1);
        
    else
        %The Bernoulli component has existence probability zero (it is
        %clutter). It will be removed by pruning
        weightB=intensity_clutter;
        
        filter_upd.tracks{Nprev_tracks+m}.eB=0;
        filter_upd.tracks{Nprev_tracks+m}.weightBLog=log(weightB);
        
        filter_upd.tracks{Nprev_tracks+m}.meanB{1}=zeros(Nx,1);
        filter_upd.tracks{Nprev_tracks+m}.covB{1}=zeros(Nx,Nx);
        filter_upd.tracks{Nprev_tracks+m}.t_ini=-1;  %We flag that these are at OOS time
        filter_upd.tracks{Nprev_tracks+m}.t_b=-1;
        filter_upd.tracks{Nprev_tracks+m}.aHis{1}=m;
        filter_upd.tracks{Nprev_tracks+m}.length=1;
        
        filter_upd.tracks{Nprev_tracks+m}.prob_length{1}=1;
        filter_upd.tracks{Nprev_tracks+m}.mean_past{1}=cell(0,1);
        filter_upd.tracks{Nprev_tracks+m}.cov_past{1}=cell(0,1);
        
        
    end
    
end


%We update all previous tracks with the measurements


for i=1:Nprev_tracks
    Nhyp_i=length(filter_retro.tracks{i}.eB);
    filter_upd.tracks{i}.t_ini=filter_retro.tracks{i}.t_ini;
    
    %We copy the birth time and length of the trajectory
    t_bi=filter_retro.tracks{i}.t_b;
    
    filter_upd.tracks{i}.t_b=t_bi;
    filter_upd.tracks{i}.length=filter_retro.tracks{i}.length;
    
    
    if(t_bi>ko)
        %Trajectory is not present at OOS time
        filter_upd.tracks{i}=filter_retro.tracks{i};
        
        for j=1:Nhyp_i
            %Misdetection hypotheses
            filter_upd.tracks{i}.weightBLog_k(j)=0; %Weight only at time step k (removing previous weight)
            %Detection hypotheses
            %We go through all measurements
            for m=1:size(z,2)
                index_hyp=j+Nhyp_i*m;
                filter_upd.tracks{i}.weightBLog_k(index_hyp)=-Inf;
            end
        end
        
    elseif(t_bi==ko)
        %Case of non-present/present trajectory
        filter_upd=retro_update_local_np(filter_retro,filter_upd,Nhyp_i,p_d,H,R,z,i,gating_threshold,Nx,Nz);
        
    else
        %The trajectory was born before ko
        %We have several cases for each possible time of death
        filter_upd=retro_update_local_before_ko(filter_retro,filter_upd,Nhyp_i,p_d,H,R,z,i,gating_threshold,Nx,Nz,ko,Lscan);
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
    
    for p=1:length(filter_retro.globHypWeight)
        %We create cost matrix for each global hypothesis.
        cost_matrix_log=-Inf(Nprev_tracks+Nnew_tracks,size(z,2));
        cost_misdetection=zeros(Nprev_tracks,1);
        
        for i=1:Nprev_tracks
            index_hyp=filter_retro.globHyp(p,i); %Hypothesis for track i in p global hypothesis
            Nhyp_i=length(filter_retro.tracks{i}.eB);
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
        
        
        
        globWeightLog_pred=log(filter_retro.globHypWeight(p));
        
        if(isempty(cost_matrix_log_trimmed))
            %One hypothesis: All targets are misdetected, no need for Murty
            opt_indices_trans=zeros(1,size(z,2));
            opt_indices_trans(fixed_indices_stay_column)=fixed_indices_stay_row;
            globWeightLog=[globWeightLog,sum(cost_misdetection)+cost_fixed+globWeightLog_pred];
            
        else
            %Number of new global hypotheses from this global hypothesis
            kbest=ceil(Nhyp_max*filter_retro.globHypWeight(p));
            
            
            %Option 1: We run Murty algorithm (making use of Hungarian algorithm to
            %solve the assginment problem). We call the function by using transpose and negative value
            %[opt_indices,nlcost]=murty(-cost_matrix_log_trimmed',kbest);
            
            %Option 2: We run MURTY algorithm (making use of Hungarian algorithm to
            %solve the assginment problem). We call the function by using transpose and negative value
            %kBest2DAssign by David F. Crouse
            %(https://github.com/USNavalResearchLaboratory/TrackerComponentLibrary)
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
            index_hyp=filter_retro.globHyp(p,i); %Hypothesis for track i in p global hypothesis
            Nhyp_i=length(filter_retro.tracks{i}.eB);
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




