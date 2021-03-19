function [X_estimate,t_b_estimate,length_estimate,tag_estimate]=PoissonMBMtarget_estimate1_tracks(filter_upd,existence_estimation_threshold,X_estimate,t_b_estimate,length_estimate,tag_estimate,k,Nx)

%Author: Angel F. Garcia Fernandez 

%We first estimate the current set of targets and their tags
%Option 1: one picks global hypothesis with highest weight and then takes the
%estimates above a certain threshold  
X_estimate_targets=[];
t_b_estimate_targets=[];
length_estimate_targets=[];
tag_estimate_targets=zeros(2,0);

if(~isempty(filter_upd.globHyp))
    globHypWeight=filter_upd.globHypWeight;
    [m,index]=max(globHypWeight);

    HypMax=filter_upd.globHyp(index,:);

    for i=1:length(HypMax)
        hyp_i=HypMax(i);
        if(hyp_i>0)
            Existence=filter_upd.tracks{i}.eB(hyp_i);
            if(Existence>existence_estimation_threshold)
       
                X_estimate_i=filter_upd.tracks{i}.meanB{hyp_i};
                tag_estimate_i=[filter_upd.tracks{i}.t_ini;filter_upd.tracks{i}.aHis{1}(1)]; %The i-th Bernoulli component has a unique tag

                X_estimate_targets=[X_estimate_targets,X_estimate_i];                        %With respect to the version that does not estimate trajectories, here we append target states as columns
                t_b_estimate_targets=[t_b_estimate_targets,filter_upd.tracks{i}.t_ini];
                length_estimate_targets=[length_estimate_targets,k-filter_upd.tracks{i}.t_ini+1];               
                tag_estimate_targets=[tag_estimate_targets,tag_estimate_i];
                
            end
        end
    end
end

%We add these estimates to the previously estimated set of trajectories
N_prev_trajectories=length(length_estimate);

%Do not forget to fill out with NaN

for i=1:size(tag_estimate_targets,2)
   %We look if there is a trajectory with the same tag 
   tag_estimate_targets_i=tag_estimate_targets(:,i);
   t_b_estimate_targets_i=t_b_estimate_targets(i);
   length_estimate_targets_i=length_estimate_targets(i);
   X_estimate_targets_i=X_estimate_targets(:,i);
   
   index_traj=find(and(tag_estimate(1,:)==tag_estimate_targets_i(1),tag_estimate(2,:)==tag_estimate_targets_i(2)));
   
   if(isempty(index_traj))
      %We have a new trajectory
      N_prev_trajectories=N_prev_trajectories+1;
      
      t_b_estimate=[t_b_estimate,t_b_estimate_targets_i];
      length_estimate=[length_estimate,length_estimate_targets_i];
      tag_estimate=[tag_estimate,tag_estimate_targets_i];
      
      X_estimate_traj=nan(Nx*length_estimate_targets_i,1);     
      X_estimate_traj(end-Nx+1:end)=X_estimate_targets_i;   
      X_estimate{N_prev_trajectories}=X_estimate_traj;
      
   else
       
       X_estimate_current=X_estimate{index_traj};
       
       X_estimate_traj=nan(Nx*length_estimate_targets_i,1);      
       X_estimate_traj(1:length(X_estimate_current))=X_estimate_current; %We fill the previous trajectory
       X_estimate_traj(end-Nx+1:end)=X_estimate_targets_i;   %We fill the current state
       
       
       X_estimate{index_traj}=X_estimate_traj;
       length_estimate(index_traj)=length_estimate_targets_i;
       
   end
    
    
end


