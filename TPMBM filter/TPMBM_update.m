function filter_upd=TPMBM_update(filter_pred,z,H,R,p_d,k,gating_threshold,intensity_clutter,Nhyp_max,Lscan)

%TPMBM update for alive trajectories

%Author: Angel F. Garcia-Fernandez
[Nz,Nx]=size(H);


Nprev_tracks=length(filter_pred.tracks);

%Poisson update
Ncom_Pois=length(filter_pred.Pois);

filter_upd.Pois=filter_pred.Pois;

%This considers that all trajectories are alive
for i=1:Ncom_Pois
    filter_upd.Pois{i}.weightPois=(1-p_d)*filter_upd.Pois{i}.weightPois;
end





%New track generation
Nnew_tracks=size(z,2); %Number of new tracks (some may have existence probability equal to zero)


for m=1:size(z,2)
    z_m=z(:,m);
    
    %We go through all Poisson components
    
    indices_new_tracks=[];
    for i=1:length(filter_pred.Pois)
        
        mean_i=filter_pred.Pois{i}.meanPois(end-Nx+1:end);
        cov_i=filter_pred.Pois{i}.covPois(end-Nx+1:end,end-Nx+1:end);
        
        S_pred_i=H*cov_i*H'+R;
        z_pred_i=H*mean_i;
        
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
        z_pred_i_max=0;
        S_pred_i_max=zeros(Nz);
        
        for i=1:length(indices_new_tracks)
            index=indices_new_tracks(i);
            mean_i=filter_pred.Pois{index}.meanPois(end-Nx+1:end);
            cov_i=filter_pred.Pois{index}.covPois(end-Nx+1:end,end-Nx+1:end);
            
            
            S_pred_i=H*cov_i*H'+R;
            z_pred_i=H*mean_i;
            
            
            weight_i=p_d*mvnpdf(z_m,z_pred_i,S_pred_i)*filter_pred.Pois{index}.weightPois;
            
            %weight_i=p_d*evaluate2DGaussianpdf(z_m,z_pred_i,S_pred_i)*filter_pred.Pois{index}.weightPois;
            
            weightB=weightB+weight_i;
            
            if(weight_i>weight_i_max)
                index_i_max=index;
                z_pred_i_max=z_pred_i;
                S_pred_i_max=S_pred_i;
                weight_i_max=weight_i;
            end
            
        end
        
        %Update of the whole trajctory (or the last Lscan time steps)
        z_pred_i=z_pred_i_max;
        S_pred_i=S_pred_i_max;
        length_k_i=filter_pred.Pois{index_i_max}.length_Pois;
        min_length=min([Lscan,length_k_i]);
        H_trajectory=[zeros(Nz,Nx*(min_length-1)),H];
        
            
        mean_i=filter_pred.Pois{index_i_max}.meanPois(end-Nx*min_length+1:end);
        cov_i=filter_pred.Pois{index_i_max}.covPois(end-Nx*min_length+1:end,end-Nx*min_length+1:end);
        
        K_pred_i=(cov_i*H_trajectory')/S_pred_i;
        P_u_i=(eye(length(mean_i))-K_pred_i*H_trajectory)*cov_i;        
        P_u_i=(P_u_i+P_u_i')/2;
        
       
        
        
        mean_u=mean_i+K_pred_i*(z_m-z_pred_i);
        
                
        meanB=filter_pred.Pois{index_i_max}.meanPois;
        meanB(end-Nx*min_length+1:end)=mean_u;
        covB=P_u_i;
        
        
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
        
        mean_j=filter_pred.tracks{i}.meanB{j};
        cov_j=filter_pred.tracks{i}.covB{j};
        eB_j=filter_pred.tracks{i}.eB(j);
        aHis_j=filter_pred.tracks{i}.aHis{j};
        
        
        %Misdetection hypotheses
        filter_upd.tracks{i}.meanB{j}=mean_j;
        filter_upd.tracks{i}.covB{j}=cov_j;
        filter_upd.tracks{i}.eB(j)=eB_j*(1-p_d)/(1-eB_j+eB_j*(1-p_d));
        filter_upd.tracks{i}.aHis{j}=[aHis_j,0];
        
        filter_upd.tracks{i}.weightBLog_k(j)=log(1-eB_j+eB_j*(1-p_d)); %Weight only at time step k (removing previous weight)
        
        
     
      
        
        %KF moments
        S_pred_j=H*cov_j(end-Nx+1:end,end-Nx+1:end)*H'+R;
        z_pred_j=H*mean_j(end-Nx+1:end);
        
        %Inverse S_pred_j
        Vs= chol(S_pred_j);
        log_det_S_pred_j= 2*log(prod(diag(Vs)));
        inv_sqrt_S= inv(Vs);
        inv_S_pred_j= inv_sqrt_S*inv_sqrt_S';
        
        
        length_k_i=filter_pred.tracks{i}.length;
        min_length=min([Lscan,length_k_i]);
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
                
                filter_upd.tracks{i}.meanB{index_hyp}=x_u;
                filter_upd.tracks{i}.covB{index_hyp}=P_u;
                filter_upd.tracks{i}.eB(index_hyp)=1;
                filter_upd.tracks{i}.aHis{index_hyp}=[aHis_j,m];
                
                %Update of the weights
                quad=-0.5*maha;
                filter_upd.tracks{i}.weightBLog_k(index_hyp)=log(eB_j*p_d)+quad-1/2*log_det_S_pred_j-Nz*log(2*pi)/2 ;
                
                
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



