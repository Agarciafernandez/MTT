function    filter_output=retro_marginal_local_before_ko(filter_retro,filter_output,i,Nhyp_i,Nx,ko,Lscan,k_is)

%This function calculates the marginal distribution of the trajectory after
%an OOS measurement.
%In this function, the trajectory was born before ko
%We have several cases for each possible time of death

t_ini=filter_retro.tracks{i}.t_ini;

if(t_ini==-1)
    %The trajectory was a first detected at the past asynchronous state, but existed from
    %before
    
    %We copy past removing last component
    
    filter_output.tracks{i}.meanB{1}=filter_retro.tracks{i}.meanB{1}(1:end-Nx);
    filter_output.tracks{i}.covB{1}=filter_retro.tracks{i}.covB{1}(1:end-Nx,1:end-Nx);
    filter_output.tracks{i}.eB(1)=filter_retro.tracks{i}.eB(1);
    filter_output.tracks{i}.aHis{1}=filter_retro.tracks{i}.aHis{1};
    filter_output.tracks{i}.prob_length{1}=filter_retro.tracks{i}.prob_length{1};
    
    %Update of the dead components
    filter_output.tracks{i}.mean_past{1}=filter_retro.tracks{i}.mean_past{1};
    filter_output.tracks{i}.cov_past{1}=filter_retro.tracks{i}.cov_past{1};
    
    %We consider as if this trajectory were born at the in sequence time
    %(it matches the required dimensions of meanB, mean_past)
    filter_output.tracks{i}.t_ini=k_is;
    
else
    filter_output.tracks{i}.t_ini=filter_retro.tracks{i}.t_ini;
    
    
    for j=1:Nhyp_i %We go through all hypotheses
        prob_length_j=filter_retro.tracks{i}.prob_length{j};
        
        if(~isempty(prob_length_j))
            %It this is varaible is empty we do not do anything. This is a local hypothesis that does not
            %belong to any global hypothesis and will be removed by pruning.
            
            
            
            
            t_f=filter_retro.tracks{i}.t_f(j);
            
            if(t_f<ko-1)
                %Trajectory is not alive at OOS time. We copy paste input to output
                
                filter_output.tracks{i}.meanB{j}=filter_retro.tracks{i}.meanB{j};
                filter_output.tracks{i}.covB{j}=filter_retro.tracks{i}.covB{j};
                filter_output.tracks{i}.eB(j)=filter_retro.tracks{i}.eB(j);
                filter_output.tracks{i}.aHis{j}=filter_retro.tracks{i}.aHis{j};
                filter_output.tracks{i}.prob_length{j}=filter_retro.tracks{i}.prob_length{j};
                
                %Update of the dead components
                filter_output.tracks{i}.mean_past{j}=filter_retro.tracks{i}.mean_past{j};
                filter_output.tracks{i}.cov_past{j}=filter_retro.tracks{i}.cov_past{j};
                
            elseif(t_f==ko-1)
                %We have a case of present/non-present
                prob_length_j2=filter_retro.tracks{i}.prob_length2{j};
                
                if(length(prob_length_j)==1 && prob_length_j==0)
                    %This corresponds to a detection hypothesis at OOS time
                    %We copy paste output, removing the last state (information is in
                    %meanB2 and covB2
                    
                    mean_j2=filter_retro.tracks{i}.meanB2{j};
                    cov_j2=filter_retro.tracks{i}.covB2{j};
                    
                    
                    %We just copy paste the information in hypothesis 2. The
                    %information in hypothesis 1 has zero weight and is removed
                    filter_output.tracks{i}.meanB{j}=mean_j2(1:end-Nx);
                    filter_output.tracks{i}.covB{j}=cov_j2(1:end-Nx,1:end-Nx);
                    filter_output.tracks{i}.prob_length{j}=prob_length_j2;
                    filter_output.tracks{i}.eB(j)=filter_retro.tracks{i}.eB(j);
                    filter_output.tracks{i}.aHis{j}=filter_retro.tracks{i}.aHis{j};
                    
                    %There are no past trajectories as this is a detection
                    filter_output.tracks{i}.mean_past{j}=cell(0,1);
                    filter_output.tracks{i}.cov_past{j}=cell(0,1);
                    
                    
                else
                    %This corresponds to a misdetection hypothesis at OOS time. We have
                    %two hypotheses for each trajectory length
                    
                    %NOTE: As this is a misdetection
                    %mean_j,cov_j and mean_j2,cov_j2 (and the past trajectories) are all equal (except at OOS time) as there has been
                    %no update. Therefore, we can just take mean_j, cov_j and merge the
                    %information in prob_length..
                    
                    mean_j=filter_retro.tracks{i}.meanB{j};
                    cov_j=filter_retro.tracks{i}.covB{j};
                    
                    
                    %Prob_length updated
                    prob_length_j_u=prob_length_j;
                    prob_length_j_u(1)=prob_length_j_u(1)+prob_length_j2; %Prob_length_j2 only contains information related to the last time step
                    
                    
                    
                    %We update the information
                    filter_output.tracks{i}.meanB{j}=mean_j(1:Nx-end);
                    filter_output.tracks{i}.covB{j}=cov_j(1:Nx-end,1:Nx-end);
                    
                    
                    filter_output.tracks{i}.prob_length{j}=prob_length_j_u;
                    filter_output.tracks{i}.eB(j)=filter_retro.tracks{i}.eB(j);
                    filter_output.tracks{i}.aHis{j}=filter_retro.tracks{i}.aHis{j};
                    
                    filter_output.tracks{i}.mean_past{j}=filter_retro.tracks{i}.mean_past{j};
                    filter_output.tracks{i}.cov_past{j}=filter_retro.tracks{i}.cov_past{j};
                    
                    
                    
                end
                
                
            elseif(t_f>=ko)
                %We have a case of present-present (though past
                %trajectories may be present/non-present)
                
                
                
                mean_j=filter_retro.tracks{i}.meanB{j};
                cov_j=filter_retro.tracks{i}.covB{j};
                prob_length_j=filter_retro.tracks{i}.prob_length{j};
                
                if(isfield(filter_retro.tracks{i},'mean_past2') && length(filter_retro.tracks{i}.prob_length2)>=j &&...
                        ~isempty(filter_retro.tracks{i}.prob_length2{j}))
                    
                    prob_length_j2=filter_retro.tracks{i}.prob_length2{j};
                else
                    prob_length_j2=zeros(size(prob_length_j));
                end
                
                
                
                
                prob_length_upd=prob_length_j+prob_length_j2;
                
                
                %We just copy paste the information in hypothesis 1. The
                %information in hypothesis 2 has zero weight (for the current time, first component in prob_length) and is removed
                filter_output.tracks{i}.meanB{j}=mean_j(1:end-Nx);
                filter_output.tracks{i}.covB{j}=cov_j(1:end-Nx,1:end-Nx);
                filter_output.tracks{i}.prob_length{j}=prob_length_upd;
                filter_output.tracks{i}.eB(j)=filter_retro.tracks{i}.eB(j);
                filter_output.tracks{i}.aHis{j}=filter_retro.tracks{i}.aHis{j};
                
                
                
                
                Npast_means=length(filter_retro.tracks{i}.mean_past{j});
                
                %Differences arise in the marginalisation of past trajectories depending if we deal with a misdetection or a detection
                %Therefore, we check this condition.
                
                if(Npast_means==0)
                    filter_output.tracks{i}.mean_past{j}=cell(0,1);
                    filter_output.tracks{i}.cov_past{j}=cell(0,1);
                    
                else
                    
                    if(filter_retro.tracks{i}.mis{j}==1)
                        %This is a misdetection hypothesis
                        
                        for p=1:Npast_means
                            t_f_p=filter_retro.tracks{i}.t_f_past{j}{p}; %We retrieve the final time step of the trajectory before OOS augmentation
                            if(t_f_p<=ko-1)
                                %We have a case of present/non-present for
                                %t_f_p=ko-1
                                %In this case we just copy paste as the past
                                %state in hypothesis 1 has not been augmented.
                                
                                %For t_f_p<ko-1, there is no state augmentation and we can just copy in the same way
                                
                                filter_output.tracks{i}.mean_past{j}{p}=filter_retro.tracks{i}.mean_past{j}{p};
                                
                                if(p<Lscan-1)
                                    filter_output.tracks{i}.cov_past{j}{p}=filter_retro.tracks{i}.cov_past{j}{p};
                                else
                                    filter_output.tracks{i}.cov_past{j}{p}=cell(0,1);
                                end
                                
                            elseif(t_f_p>=ko)
                                %We have a case of present/present. We need to
                                %remove the last state as it has been augmented
                                
                                filter_output.tracks{i}.mean_past{j}{p}=filter_retro.tracks{i}.mean_past{j}{p}(1:end-Nx);
                                
                                if(p<Lscan-1)
                                    filter_output.tracks{i}.cov_past{j}{p}=filter_retro.tracks{i}.cov_past{j}{p}(1:end-Nx,1:end-Nx);
                                else
                                    filter_output.tracks{i}.cov_past{j}{p}=cell(0,1);
                                end
                            end
                        end
                        
                    else
                        %This a detection hypothesis
                        
                        for p=1:Npast_means
                            t_f_p=filter_retro.tracks{i}.t_f_past{j}{p}; %We retrieve the final time step of the trajectory before OOS augmentation
                            if(t_f_p==ko-1)
                                %We have a case of present/non-present for
                                %t_f_p=ko-1
                                %The detection mean and covariance are in
                                %hypothesis 2
                                
                                
                                filter_output.tracks{i}.mean_past{j}{p}=filter_retro.tracks{i}.mean_past2{j}{p}(1:end-Nx);
                                filter_output.tracks{i}.cov_past{j}{p}=filter_retro.tracks{i}.cov_past2{j}{p}(1:end-Nx,1:end-Nx);
                                
                            elseif(t_f_p>=ko)
                                %We have a case of present/present. We need to
                                %remove the last state as it has been augmented
                                %Stored in hypothesis one
                                
                                filter_output.tracks{i}.mean_past{j}{p}=filter_retro.tracks{i}.mean_past{j}{p}(1:end-Nx);
                                filter_output.tracks{i}.cov_past{j}{p}=filter_retro.tracks{i}.cov_past{j}{p}(1:end-Nx,1:end-Nx);
                                
                            else
                                % Case t_f_p<ko-1, there is no state
                                % augmentation and we can just copy the past
                                % mean and covariance
                                filter_output.tracks{i}.mean_past{j}{p}=filter_retro.tracks{i}.mean_past{j}{p};
                                filter_output.tracks{i}.cov_past{j}{p}=filter_retro.tracks{i}.cov_past{j}{p};
                                
                                
                                
                            end
                        end
                    end
                    
                end
                
                
            end
        end
    end
    
end
