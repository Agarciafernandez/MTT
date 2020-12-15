function filter_upd=TMBM_all_update(filter_pred,z,H,R,p_d,k,gating_threshold,intensity_clutter,Nhyp_max,Lscan,T_alive)

%TMBM update for alive trajectories. Same as TPMBM but with Poisson
%intensity equal to zero (with an improvement for this the multi-Bernoulli case done by Yuxuan Xia for the MBM filter for sets of targets)

%Author: Angel F. Garcia-Fernandez

[Nz,Nx]=size(H);

Nprev_tracks=length(filter_pred.tracks);


%We update all previous tracks with the measurements

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
        
        prob_length_j=filter_pred.tracks{i}.prob_length{j};
        prob_alive_j=prob_length_j(1);
        
        if(prob_alive_j>T_alive)
            %Trajectory is considered alive
            
            %We compute some required terms
            cross_p_d=p_d*prob_alive_j;
            
            
            
            %Misdetection hypotheses
            filter_upd.tracks{i}.meanB{j}=mean_j;
            filter_upd.tracks{i}.covB{j}=cov_j;
            filter_upd.tracks{i}.eB(j)=eB_j*(1-cross_p_d)/(1-eB_j+eB_j*(1-cross_p_d));
            filter_upd.tracks{i}.aHis{j}=[aHis_j,0];
            filter_upd.tracks{i}.weightBLog_k(j)=log(1-eB_j+eB_j*(1-cross_p_d)); %Weight only at time step k (removing previous weight)
            
            %Update the length probabilities
            prob_length_j_upd=prob_length_j;
            prob_length_j_upd(1)=prob_length_j_upd(1)*(1-p_d);
            prob_length_j_upd=prob_length_j_upd/(1-cross_p_d);
            filter_upd.tracks{i}.prob_length{j}=prob_length_j_upd;
            
            %Update the dead components            
            filter_upd.tracks{i}.mean_past{j}=filter_pred.tracks{i}.mean_past{j};
            
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
                    p_z_log=quad-1/2*log_det_S_pred_j-Nz*log(2*pi)/2;
                    filter_upd.tracks{i}.weightBLog_k(index_hyp)=log(eB_j*cross_p_d)+p_z_log; %We use cross_p_d
                    
                    %Update the length probabilities and mean_past
                    filter_upd.tracks{i}.prob_length{index_hyp}=1;
                    filter_upd.tracks{i}.mean_past{index_hyp}=cell(0,1);       
                    
                else
                    filter_upd.tracks{i}.weightBLog_k(index_hyp)=-Inf;
                    
                    
                end
                
            end
        else
            %Trajectory is considered dead
            %Misdetection hypotheses
            filter_upd.tracks{i}.meanB{j}=mean_j;
            filter_upd.tracks{i}.covB{j}=cov_j;
            filter_upd.tracks{i}.eB(j)=eB_j;
            filter_upd.tracks{i}.aHis{j}=[aHis_j,0];
            filter_upd.tracks{i}.weightBLog_k(j)=0; %Weight only at time step k (removing previous weight)
            filter_upd.tracks{i}.prob_length{j}=prob_length_j;
            
            %Update of the dead components           
            filter_upd.tracks{i}.mean_past{j}=filter_pred.tracks{i}.mean_past{j};
            
            
            %Detection hypotheses
            %We go through all measurements
            for m=1:size(z,2)
                index_hyp=j+Nhyp_i*m;
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
    filter_upd.globHyp=ones(1,0);
    
else
    
    globWeightLog=[];
    globHyp=[];
    
    for p=1:length(filter_pred.globHypWeight)
        %We create cost matrix for each global hypothesis.
        cost_matrix_log=-Inf(Nprev_tracks,size(z,2));
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
                
                %We remove the weight of the clutter intensity to use Murty (Later this weight is added).
                cost_matrix_log(i,1:length(indices_c))=weights_log-log(intensity_clutter);

                cost_misdetection(i)=filter_upd.tracks{i}.weightBLog_k(index_hyp);
                
                
            end
        end
        
       %Clutter hypotheses 
        cost_temp = -Inf(Nprev_tracks);
        cost_temp(logical(eye(Nprev_tracks))) = cost_misdetection;
        
        %Further process cost_matrix_log, remove rows and columns that are
        %all -inf
        indices_not_inf = ~isinf(cost_matrix_log);
        indices_valid_row_boolean = sum(indices_not_inf,2) > 0;
        
        indices_valid_column_boolean = sum(indices_not_inf,1) > 0;
        indices_valid_column = find(indices_valid_column_boolean);
        num_valid_meas = length(indices_valid_column);
        
        cost_matrix_log = [cost_matrix_log(indices_valid_row_boolean,indices_valid_column_boolean) cost_temp(indices_valid_row_boolean,indices_valid_row_boolean)];
%         cost_matrix_log = [cost_matrix_log cost_temp];

        globWeightLog_pred=log(filter_pred.globHypWeight(p));
        
        if(isempty(cost_matrix_log))
            %One hypothesis: All targets are misdetected, no need for Murty                 
            opt_indices=zeros(1,Nprev_tracks);
            globWeightLog=[globWeightLog,sum(cost_misdetection)+globWeightLog_pred+log(intensity_clutter)*size(z,2)];
      
        else
        
        
        %Number of new global hypotheses from this global hypothesis
        kbest=ceil(Nhyp_max*filter_pred.globHypWeight(p));
        
        %We run Murty algorithm (making use of Hungarian algorithm to
        %solve the assginment problem). We call the function by using transpose and negative value
        [opt_indices_temp,nlcost]= murty(-cost_matrix_log,kbest);
        
        
            
            %Option 2: We run MURTY algorithm (making use of Hungarian algorithm to
            %solve the assginment problem). We call the function by using transpose and negative value
            %kBest2DAssign by David F. Crouse
            [opt_indices_temp,~,nlcost]=kBest2DAssign(-cost_matrix_log,kbest);
            opt_indices_temp=opt_indices_temp';
            nlcost=nlcost';
        
        opt_indices_temp(opt_indices_temp>num_valid_meas) = 0;
        
        opt_indices = zeros(size(opt_indices_temp,1),Nprev_tracks);
        opt_indices(:,indices_valid_row_boolean) = opt_indices_temp;
        
        for i = 1:size(opt_indices,1)
            for j = 1:Nprev_tracks
                if opt_indices(i,j)~=0
                    opt_indices(i,j) = indices_valid_column(opt_indices(i,j));
%                     temp = indices_valid_column(opt_indices(i,j));
%                     opt_indices(i,j) = opt_indices(i,j) + temp - sum(indices_valid_column_boolean(1:temp));
                end
            end
        end
        
        
        globWeightLog=[globWeightLog,-nlcost+globWeightLog_pred+log(intensity_clutter)*size(z,2)+sum(cost_misdetection(~indices_valid_row_boolean))];
        
        end
        
        %We write the corresponding global hypothesis
        globHypProv=zeros(size(opt_indices,1),Nprev_tracks);
        
        for i=1:Nprev_tracks
            index_hyp=filter_pred.globHyp(p,i); %Hypothesis for track i in p global hypothesis
            Nhyp_i=length(filter_pred.tracks{i}.eB);
            for j=1:size(opt_indices,1)
                if opt_indices(j,i) == 0
                    index=index_hyp;
                else
                    index=index_hyp+Nhyp_i*opt_indices(j,i);
                end
                globHypProv(j,i)=index;
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





