function filter_output=retro_marginal_local_np(filter_retro,filter_output,i,Nhyp_i,Nx)

%This function calculates the marginal distribution of the trajectory after
%an OOS measurement. The hypothesis is non-present/present.

for j=1:Nhyp_i %We go through all hypotheses
    
    prob_length_j=filter_retro.tracks{i}.prob_length{j};
    
    if(isempty(prob_length_j))
        %Empty 
        %We do not do anything. This is a local hypothesis that does not
        %belong to any global hypothesis and will be removed by pruning.

    
    elseif(length(prob_length_j)==1 && prob_length_j==0)
        %This corresponds to a detection hypothesis at OOS time
        %We copy paste output, removing the last state (information is in
        %meanB2 and covB2
        mean_j2=filter_retro.tracks{i}.meanB2{j};
        cov_j2=filter_retro.tracks{i}.covB2{j};
        prob_length_j2=filter_retro.tracks{i}.prob_length2{j};
        
        
        %We just copy paste the information in hypothesis 2. The
        %information in hypothesis 1 has zero weight and is removed
        filter_output.tracks{i}.meanB{j}=mean_j2(1:end-Nx);
        filter_output.tracks{i}.covB{j}=cov_j2(1:end-Nx,1:end-Nx);
        filter_output.tracks{i}.prob_length{j}=prob_length_j2;
        filter_output.tracks{i}.eB(j)=filter_retro.tracks{i}.eB(j);
        filter_output.tracks{i}.aHis{j}=filter_retro.tracks{i}.aHis{j};
        filter_output.tracks{i}.t_ini=filter_retro.tracks{i}.t_ini;
        
        %Now we consider past trajectories
        if(isfield(filter_retro.tracks{i},'mean_past2') && length(filter_retro.tracks{i}.mean_past2)>=j )
            Npast_means=length(filter_retro.tracks{i}.mean_past2{j});
            for p=1:Npast_means
                filter_output.tracks{i}.mean_past{j}{p}=filter_retro.tracks{i}.mean_past2{j}{p}(1:end-Nx);
                filter_output.tracks{i}.cov_past{j}{p}=filter_retro.tracks{i}.cov_past2{j}{p}(1:end-Nx,1:end-Nx);
            end
        else
            filter_output.tracks{i}.mean_past{j}=cell(0,1);
            filter_output.tracks{i}.cov_past{j}=cell(0,1);
        end
        
    else
        %This corresponds to a misdetection hypothesis at OOS time. We have
        %two hypotheses for each trajectory length
        
        %NOTE: As this is a misdetection 
        %mean_j,cov_j and mean_j2,cov_j2 (and the past trajectories) are all equal (except at OOS time) as there has been
        %no update. Therefore, we can just take mean_j, cov_j and merge the
        %information in prob_length.
        
        mean_j=filter_retro.tracks{i}.meanB{j};
        cov_j=filter_retro.tracks{i}.covB{j};
        
        prob_length_j2=filter_retro.tracks{i}.prob_length2{j};

        %Prob_length updated
        prob_length_j_u=(prob_length_j+prob_length_j2);
        
          %We update the information
        filter_output.tracks{i}.meanB{j}=mean_j;
        filter_output.tracks{i}.covB{j}=cov_j;
        filter_output.tracks{i}.prob_length{j}=prob_length_j_u;
        filter_output.tracks{i}.eB(j)=filter_retro.tracks{i}.eB(j);
        filter_output.tracks{i}.aHis{j}=filter_retro.tracks{i}.aHis{j};
        
        filter_output.tracks{i}.mean_past{j}=filter_retro.tracks{i}.mean_past{j};
        filter_output.tracks{i}.cov_past{j}=filter_retro.tracks{i}.cov_past{j};
                filter_output.tracks{i}.t_ini=filter_retro.tracks{i}.t_ini;

        
        %OLD
        
%         mean_j=filter_retro.tracks{i}.meanB{j};
%         cov_j=filter_retro.tracks{i}.covB{j};
%         
%         
%         mean_j2=filter_retro.tracks{i}.meanB2{j}(1:end-Nx);
%         cov_j2=filter_retro.tracks{i}.covB2{j}(1:end-Nx,1:end-Nx);
%         
%         prob_length_j2=filter_retro.tracks{i}.prob_length2{j};
%         
%         %Prob_length updated
%         prob_length_j_u=(prob_length_j+prob_length_j2);
%                 
%         [mean_j_u,cov_j_u]=fitGaussian(prob_length_j(1),prob_length_j2(1),mean_j,cov_j,mean_j2,cov_j2);
%         
%         %We update the informatio
%         filter_output.tracks{i}.meanB{j}=mean_j_u;
%         filter_output.tracks{i}.covB{j}=cov_j_u;
%         filter_output.tracks{i}.prob_length{j}=prob_length_j_u;
%         filter_output.tracks{i}.eB(j)=filter_retro.tracks{i}.eB(j);
%         filter_output.tracks{i}.aHis{j}=filter_retro.tracks{i}.aHis{j};
%         
%         %We update past means
%         
%         Npast_means=length(filter_retro.tracks{i}.mean_past{j});
%         for p=1:Npast_means
%             mean_j=filter_retro.tracks{i}.mean_past{j}{p};
%             cov_j=filter_retro.tracks{i}.cov_past{j}{p};
%             
%             mean_j2=filter_retro.tracks{i}.mean_past2{j}{p}(1:end-Nx);
%             cov_j2=filter_retro.tracks{i}.cov_past2{j}{p}(1:end-Nx,1:end-Nx);
%             
%             
%             [mean_j_u,cov_j_u]=fitGaussian(prob_length_j(p+1),prob_length_j2(p+1),mean_j,cov_j,mean_j2,cov_j2);
%             
%             filter_output.tracks{i}.mean_past{j}{p}=mean_j_u;
%             filter_output.tracks{i}.cov_past{j}{p}=cov_j_u;
%             
%         end
%         
        
    end
    
    
    
end
end


% function [mean_j_u,cov_j_u]=fitGaussian(p1u,p2u,mean_j,cov_j,mean_j2,cov_j2)
% %Moment matching for two Gaussians
% ptot=p1u+p2u;
% p1=p1u/ptot;
% p2=p2u/ptot;
% 
% %We perform moment matching to obtain the update
% mean_j_u=p1*mean_j+p2*mean_j2;
% cov_j_u=p1*(cov_j+(mean_j-mean_j_u)*(mean_j-mean_j_u)')+p2*(cov_j2+(mean_j2-mean_j_u)*(mean_j2-mean_j_u)');
% 
% end