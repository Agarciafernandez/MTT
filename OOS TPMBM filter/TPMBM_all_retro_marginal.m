function filter_output=TPMBM_all_retro_marginal(filter_upd,Nx,Lscan,k_is)

ko=filter_upd.ko;

%We copy the information on global hypotheses from input to output

filter_output.globHypWeight=filter_upd.globHypWeight;
filter_output.globHyp=filter_upd.globHyp;



%This function performs the marginalisation on the TPMBM density after the
%update step with the OOS measurement


%We start with the PPP
Ncom_Pois=length(filter_upd.Pois);
for i=1:Ncom_Pois
    
    filter_output.Pois{i}.t_bPois=filter_upd.Pois{i}.t_bPois; %To mark that it is an OOS state
    filter_output.Pois{i}.length_Pois=filter_upd.Pois{i}.length_Pois;
    
    if(filter_upd.Pois{i}.is_oos)
        mean_i=filter_upd.Pois{i}.meanPois;
        if(length(mean_i)>Nx)
            cov_i=filter_upd.Pois{i}.covPois;
            %We remove the last state in the mean and covariance matrix
            filter_output.Pois{i}.meanPois=mean_i(1:end-Nx);
            filter_output.Pois{i}.covPois=cov_i(1:end-Nx,1:end-Nx);
            filter_output.Pois{i}.weightPois=filter_upd.Pois{i}.weightPois;
            
            
        else
            %The component only has an OOS state (so this component is
            %removed) (It will be effectively removed after pruning
            filter_output.Pois{i}.meanPois=cell(0,1);
            filter_output.Pois{i}.covPois=cell(0,1);
            filter_output.Pois{i}.weightPois=0;
            
        end
        
    else
        %We leave this component unchanged
        filter_output.Pois{i}.meanPois=filter_upd.Pois{i}.meanPois;
        filter_output.Pois{i}.covPois=filter_upd.Pois{i}.covPois;
        filter_output.Pois{i}.weightPois=filter_upd.Pois{i}.weightPois;
        
        
    end
end



Ntracks=length(filter_upd.tracks);

if(Ntracks>0)
    
    for i=1:Ntracks
        Nhyp_i=length(filter_upd.tracks{i}.eB);
        
        
        %t_ini,t_b and length do not change
        t_bi=filter_upd.tracks{i}.t_b;
        
        filter_output.tracks{i}.t_b=t_bi;
        filter_output.tracks{i}.length=filter_upd.tracks{i}.length;
        
        if(t_bi>ko)
            %Trajectory is not present at OOS time so it is copied
            %unchanged
            filter_output.tracks{i}=filter_upd.tracks{i};
            
        elseif(t_bi==ko)
            %Scenario of non-present/present
            filter_output=retro_marginal_local_np(filter_upd,filter_output,i,Nhyp_i,Nx);
        elseif(t_bi==-1)
            %Trajectory was created at OOS time. It does not exist at other
            %times.
            %We set the existence probability equal to zero, so that it is pruned
            
            filter_output.tracks{i}=filter_upd.tracks{i};
            filter_output.tracks{i}.eB=0;
            
            
        else
            %The trajectory was born before ko
            %We have several cases for each possible time of death
            filter_output=retro_marginal_local_before_ko(filter_upd,filter_output,i,Nhyp_i,Nx,ko,Lscan,k_is);

        end
        
    end
end