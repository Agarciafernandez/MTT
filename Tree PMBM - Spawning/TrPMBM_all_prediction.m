function filter_pred=TrPMBM_all_prediction(filter_upd,F,Q,p_s,weights_b,means_b,covs_b,Lscan,k,T_alive,F1,F2,p_s1,p_s2,d_spawning)
%Author: Angel F. Garcia-Fernandez 

%This function works for two spawning modes, see Scenario_spawning. It is
%straightforward to add more spawning modes.

%Prediction for Poisson component (we only keep track of alive trajectories
%in the Poisson component)
Ncom=length(filter_upd.Pois);
Nx=size(F,1);

for i=1:Ncom
    filter_pred.Pois{i}.weightPois=p_s*filter_upd.Pois{i}.weightPois;
    filter_pred.Pois{i}.meanPois=[filter_upd.Pois{i}.meanPois;F*filter_upd.Pois{i}.meanPois(end-Nx+1:end)];
    
    %Prediction of the covariance matrix
    length_k_i=filter_upd.Pois{i}.length_Pois;
    if(Lscan>1)
        min_length=min([Lscan-1,length_k_i]);
        cov_i=filter_upd.Pois{i}.covPois(end-Nx*min_length+1:end,end-Nx*min_length+1:end);
        cov_i_pred=F*cov_i(end-Nx+1:end,end-Nx+1:end)*F'+Q;
        F_mode=zeros(Nx,size(cov_i,1));
        F_mode(end-3:end,end-3:end)=F;
        cross=F_mode*cov_i;
        cov_upd_i=[cov_i,cross';cross,cov_i_pred];
    else
        cov_i=filter_upd.Pois{i}.covPois;
        cov_upd_i=F*cov_i*F'+Q;
        
    end
    
    filter_pred.Pois{i}.covPois=cov_upd_i;
    %Now we predict the start time and length of the trajectory
    filter_pred.Pois{i}.t_bPois=filter_upd.Pois{i}.t_bPois;
    filter_pred.Pois{i}.length_Pois=filter_upd.Pois{i}.length_Pois+1;
end
%We add the PHD of new born trajectories
for i=1:length(weights_b)
    filter_pred.Pois{Ncom+i}.weightPois=weights_b(i);
    filter_pred.Pois{Ncom+i}.meanPois=means_b(:,i);
    filter_pred.Pois{Ncom+i}.covPois=squeeze(covs_b(:,:,i));
    filter_pred.Pois{Ncom+i}.t_bPois=k+1; %This birth component is born at time k+1
    filter_pred.Pois{Ncom+i}.length_Pois=1;
end




%Prediction for Bernoulli components
globHyp_upd=filter_upd.globHyp;

filter_pred.globHypWeight=filter_upd.globHypWeight;




Ntracks=length(filter_upd.tracks);
index_tracks_alive=false(Ntracks,1);

if(Ntracks>0)
    for i=1:Ntracks
        Nhyp_i=length(filter_upd.tracks{i}.eB);
        filter_pred.tracks{i}.t_ini=filter_upd.tracks{i}.t_ini;
        %We predict the birth time and length of the trajectory
        filter_pred.tracks{i}.t_b=filter_upd.tracks{i}.t_b;
        filter_pred.tracks{i}.length=filter_upd.tracks{i}.length+1;
        
        genealogy_i=filter_upd.tracks{i}.genealogy;
        
        filter_pred.tracks{i}.genealogy=[genealogy_i;1];
        
        for j=1:Nhyp_i %We go through all hypotheses
            
            
            %We first check if the trajectory is alive (can be alive)
            prob_length_j=filter_upd.tracks{i}.prob_length{j};
            
            if(prob_length_j(1)>T_alive)
                
                index_tracks_alive(i)=1;
                
                prob_alive=p_s*prob_length_j(1);
                
                %Prediction of the mean
                filter_pred.tracks{i}.meanB{j}=[filter_upd.tracks{i}.meanB{j};F*filter_upd.tracks{i}.meanB{j}(end-Nx+1:end)];
                
                %Prediction of the covariance matrix
                length_k_i=filter_upd.tracks{i}.length;
                
                if(Lscan>1)
                    min_length=min([Lscan-1,length_k_i]);
                    cov_i=filter_upd.tracks{i}.covB{j}(end-Nx*min_length+1:end,end-Nx*min_length+1:end);
                    cov_i_pred=F*cov_i(end-Nx+1:end,end-Nx+1:end)*F'+Q;
                    F_mode=zeros(Nx,size(cov_i,1));
                    F_mode(end-3:end,end-3:end)=F;
                    cross=F_mode*cov_i;
                    cov_upd_i=[cov_i,cross';cross,cov_i_pred];
                else
                    cov_i=filter_upd.tracks{i}.covB{j};
                    cov_upd_i=F*cov_i*F'+Q;
                end
                
                
                filter_pred.tracks{i}.covB{j}=cov_upd_i;
                filter_pred.tracks{i}.eB(j)=filter_upd.tracks{i}.eB(j); %Difference w.r.t. alive trajectories
                filter_pred.tracks{i}.aHis{j}=filter_upd.tracks{i}.aHis{j};
                
                %We update the probability of length
                filter_pred.tracks{i}.prob_length{j}=[prob_alive;(1-p_s)*prob_length_j(1);filter_upd.tracks{i}.prob_length{j}(2:end)];
                
                
                %We update the past means of the trajectories
                filter_pred.tracks{i}.mean_past{j}{1}=filter_upd.tracks{i}.meanB{j};
                %Npast_means=length(prob_length_j)-1;
                Npast_means=length(filter_upd.tracks{i}.mean_past{j});
                for p=1:Npast_means
                    filter_pred.tracks{i}.mean_past{j}{p+1}=filter_upd.tracks{i}.mean_past{j}{p};
                end
                
                
                
                
            else
                %We just copy paste the previous Bernoulli component
                
                filter_pred.tracks{i}.meanB{j}=filter_upd.tracks{i}.meanB{j};
                filter_pred.tracks{i}.covB{j}=filter_upd.tracks{i}.covB{j};
                filter_pred.tracks{i}.eB(j)=filter_upd.tracks{i}.eB(j);
                filter_pred.tracks{i}.aHis{j}=filter_upd.tracks{i}.aHis{j};
                filter_pred.tracks{i}.prob_length{j}=prob_length_j;
                filter_pred.tracks{i}.mean_past{j}=filter_upd.tracks{i}.mean_past{j};
            end
            
            
        end
        
    end
else
    filter_pred.tracks=cell(0,1);
    filter_pred.globHyp=[];
    filter_pred.globHypWeight=[];
end


%Now we consider the newly spawning Bernoullis
Ntracks_alive=sum(index_tracks_alive);

%We repeat the global hypothesis matrix twice for the alive tracks
%(branches), one for each spawning mode
filter_pred.globHyp=[globHyp_upd,globHyp_upd(:,index_tracks_alive),globHyp_upd(:,index_tracks_alive)];

%We predict tree_index variable
filter_pred.N_tree=filter_upd.N_tree;

tree_index_upd=filter_upd.tree_index;

filter_pred.tree_index=[tree_index_upd;...
    tree_index_upd(index_tracks_alive);...
    tree_index_upd(index_tracks_alive)];

index_output=1; %Index to write the output

for i=1:Ntracks
    
    Nhyp_i=length(filter_upd.tracks{i}.eB);

    
    if(index_tracks_alive(i))
        
        genealogy_i=filter_upd.tracks{i}.genealogy;
        index_sp1=index_output+Ntracks;
        index_sp2=index_output+Ntracks+Ntracks_alive;
        
        %Metadata spawning mode 1
        filter_pred.tracks{index_sp1}.genealogy=[genealogy_i;2];
        filter_pred.tracks{index_sp1}.t_ini=k+1;
        filter_pred.tracks{index_sp1}.t_b=k+1;
        filter_pred.tracks{index_sp1}.length=1;
        
        %Metadata spawning mode 2
        filter_pred.tracks{index_sp2}.genealogy=[genealogy_i;3];
        filter_pred.tracks{index_sp2}.t_ini=k+1;
        filter_pred.tracks{index_sp2}.t_b=k+1;
        filter_pred.tracks{index_sp2}.length=1;
        
        for j=1:Nhyp_i %We go through all hypotheses
            %We need to calculate probability of existence and single
            %target hypothesis
            prob_length_j=filter_upd.tracks{i}.prob_length{j};
            eBj=filter_upd.tracks{i}.eB(j);
            
            mean_j_act=filter_upd.tracks{i}.meanB{j}(end-Nx+1:end); %Mean of the state at the current time
            cov_j_act=filter_upd.tracks{i}.covB{j}(end-Nx+1:end,end-Nx+1:end); %Covariance of the state at the current time
            
            
            vel=sqrt(mean_j_act(2)^2+mean_j_act(4)^2);
            unit_vector=[mean_j_act(2)/vel,0,mean_j_act(4)/vel,0]';
            
            d_offset1=d_spawning*[-unit_vector(3),0,unit_vector(1),0]';
            d_offset2=-d_offset1; %Signed changed w.r.t. spawning mode 1

            
            %Spawning mode 1
            eB1=p_s1*prob_length_j(1)*eBj;
            mean1=d_offset1+F1*mean_j_act;
            cov1=F1*cov_j_act*F1'+Q;
            
            
            filter_pred.tracks{index_sp1}.meanB{j}=mean1;
            filter_pred.tracks{index_sp1}.covB{j}=cov1;
            filter_pred.tracks{index_sp1}.eB(j)=eB1;
            filter_pred.tracks{index_sp1}.aHis{j}=filter_upd.tracks{i}.aHis{j};
            filter_pred.tracks{index_sp1}.prob_length{j}=1;
            filter_pred.tracks{index_sp1}.mean_past{j}=cell(0,1);
            
            
            %Spawning mode 2
            eB2=p_s2*prob_length_j(1)*eBj;
            mean2=d_offset2+F2*mean_j_act;
            cov2=F2*cov_j_act*F2'+Q;

            filter_pred.tracks{index_sp2}.meanB{j}=mean2;
            filter_pred.tracks{index_sp2}.covB{j}=cov2;
            filter_pred.tracks{index_sp2}.eB(j)=eB2;
            filter_pred.tracks{index_sp2}.aHis{j}=filter_upd.tracks{i}.aHis{j};
            filter_pred.tracks{index_sp2}.prob_length{j}=1;
            filter_pred.tracks{index_sp2}.mean_past{j}=cell(0,1);
                       
            
        end
        
        index_output=index_output+1;
        
    end
end


