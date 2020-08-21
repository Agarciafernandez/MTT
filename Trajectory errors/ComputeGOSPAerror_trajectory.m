function    [squared_gospa,gospa_loc,gospa_mis,gospa_fal]=ComputeGOSPAerror_trajectory(X_estimate,t_b_estimate, length_estimate,X_truth,t_birth,t_death,c_gospa,k_end,Nx)


%This function computes the (squared) error of the alive trajectories and its
%decomposition using the sum of the GOSPA metric at each time step (normalised by the current time window, which is from 1
%to k_end)
%Author: Angel F. Garcia-Fernandez


squared_gospa=0;
gospa_loc=0;
gospa_mis=0;
gospa_fal=0;

for k=1:k_end
    %We obtain X (we only consider alive trajectories at the current time k_end)
    alive_targets_index=and(and(and(k>=t_birth,k<t_death),k_end>=t_birth),k_end<t_death);
    Nalive_targets=sum(alive_targets_index);
    alive_targets_pos_index=[alive_targets_index;false(1,length(t_birth));alive_targets_index;false(1,length(t_birth))];
    X_k=reshape(X_truth(alive_targets_pos_index(:),k),2,Nalive_targets);

    %We obtain Y
    alive_targets_index=find(and(k>=t_b_estimate,k<=t_b_estimate+length_estimate-1));
    Y_k=zeros(Nx,length(alive_targets_index));   
    for i=1:length(alive_targets_index)    
        index_loc=k-t_b_estimate(alive_targets_index(i))+1;       
        Y_k(:,i)=X_estimate{alive_targets_index(i)}((index_loc-1)*Nx+1:index_loc*Nx);      
    end
    %Rremove nan values in the estimate (which are created by holes in the
    %trajectories)
    index_nonan=~isnan(Y_k(1,:)); 
    Y_k_pos=Y_k([1 3],index_nonan);
    
    
    [d_gospa, ~, decomp_cost] = GOSPA(X_k, Y_k_pos, 2, c_gospa, 2); 
    squared_gospa=squared_gospa+d_gospa^2;
    gospa_loc=gospa_loc+decomp_cost.localisation;
    gospa_mis=gospa_mis+decomp_cost.missed;
    gospa_fal=gospa_fal+decomp_cost.false;
end

squared_gospa=squared_gospa/k_end;
gospa_loc=gospa_loc/k_end;
gospa_mis=gospa_mis/k_end;
gospa_fal=gospa_fal/k_end;
