function filter_retro=TPMBM_all_retro(filter_upd,to,tk_is,k,k_is,q,mean_a_p,mean_a_v,P_a_pp,P_a_vv,P_a_pv,lambda,mu,Lscan)

%This function performs the retrodiction step for the TPMBM filter to
%process an OOS measurement at time to

%Many of the parameters in the retrodicted PMBM density and in the posterior PMBM are
%alike, so we just modify the requred ones.
filter_retro=filter_upd;

Nx=4;


ko=find(tk_is>to,1);

delta_t1=to-tk_is(ko-1); %Time to the previous in-sequence measurement
delta_t2=tk_is(ko)-to; %Time to the next in-sequence measurement

delta_t_total=delta_t1+delta_t2;

%Probabilities of survival within each interval (not to confuse with the
%probability of survival with in-sequence measurements)
p_s1=(exp(-mu*delta_t1)-exp(-mu*(delta_t_total)))/(1-exp(-mu*delta_t_total));
p_s2=(exp(-mu*delta_t2)-exp(-mu*(delta_t_total)))/(1-exp(-mu*delta_t_total));


%Transition matrices in each interval
F1=kron(eye(2),[1 delta_t1;0 1]);
Q1=q*kron(eye(2),[delta_t1^3/3 delta_t1^2/2; delta_t1^2/2 delta_t1]);
F2=kron(eye(2),[1 delta_t2;0 1]);
Q2=q*kron(eye(2),[delta_t2^3/3 delta_t2^2/2; delta_t2^2/2 delta_t2]);

%Birth model at OOS time
[mean_b_1,cov_b_1]=BirthParameters_cd(q,delta_t1,mean_a_p,mean_a_v,P_a_pp,P_a_vv,P_a_pv,mu);
weights_b=lambda/mu*(1-exp(-mu*delta_t1))*(1-exp(-mu*delta_t2));


%We precompute some more matrices
%Matrices for hypotheses present-present
Kpp=Q1*F2'/(F2*Q1*F2'+Q2);
Qpp=Q1-Kpp*F2*Q1;
Fpp_prov=F1-Kpp*F2*F1;

%Matrices for hypotheses nonpresent-present
Knp=cov_b_1*F2'/(F2*cov_b_1*F2'+Q2);
Qnp=cov_b_1-Knp*F2*cov_b_1;
xb_prov_np=(eye(Nx)-Knp*F2)*mean_b_1;

%Retrodiction for the PPP (which only considers alive trajectories)
Ncom=length(filter_retro.Pois);
index_extra=1; %Index to store newly generated components

%We store p_s1,p_s2 and ko parameters separatedly (as all newly create mixtures have
%two components with these weights)
filter_retro.p_s1=p_s1;
filter_retro.p_s2=p_s2;
filter_retro.ko=ko;

for i=1:Ncom
    
    
    t_bi=filter_retro.Pois{i}.t_bPois;
    
    %We recall that the OOS is within  the Lscan window
    if(t_bi<ko)
        %Trajectory is alive at both time steps
        
        mean_i=filter_retro.Pois{i}.meanPois;
        cov_i=filter_retro.Pois{i}.covPois;
        
        [mean_o,cov_o]=retro_pp(mean_i,cov_i, Kpp,Qpp,Fpp_prov,Nx,ko,k_is);
        
        
        %We write the parameters in the output
        filter_retro.Pois{i}.meanPois=mean_o;
        filter_retro.Pois{i}.covPois=cov_o;
        %We add the flag that this component has an OOS state
        filter_retro.Pois{i}.is_oos=1;
        
        
    elseif(t_bi==ko)
        %Trajectory is nonpresent-present at the time interval
        %We have two hypotheses
        %The one that does not change is the same as before, so we just
        %have (with modified weight
        filter_retro.Pois{i}.is_oos=0;
        filter_retro.Pois{i}.weightPois=(1-p_s2)*filter_upd.Pois{i}.weightPois;
        
        %We compute the newly added component
        
        
        mean_i=filter_upd.Pois{i}.meanPois;
        cov_i=filter_upd.Pois{i}.covPois;
        
        [mean_o,cov_o]=retro_np(mean_i,cov_i,Knp,Qnp,xb_prov_np,Nx);
        
        
        %We write the parameters in the output
        filter_retro.Pois{Ncom+index_extra}.meanPois=mean_o;
        filter_retro.Pois{Ncom+index_extra}.covPois=cov_o;
        %We add the flag that this component has an OOS state
        filter_retro.Pois{Ncom+index_extra}.is_oos=1;
        filter_retro.Pois{Ncom+index_extra}.weightPois=p_s2*filter_upd.Pois{i}.weightPois;
        filter_retro.Pois{Ncom+index_extra}.t_bPois=-1; %To mark that it is an OOS state
        filter_retro.Pois{Ncom+index_extra}.length_Pois=1; %To mark that it is an OOS state
        
        
        index_extra=index_extra+1;
        
    else
        %We add the flag that this component does not have an OOS state
        filter_retro.Pois{i}.is_oos=0;
        
    end
    
end

%We add the PHD of trajectories born in the OOS time
filter_retro.Pois{Ncom+index_extra}.weightPois=weights_b;
filter_retro.Pois{Ncom+index_extra}.meanPois=mean_b_1;
filter_retro.Pois{Ncom+index_extra}.covPois=cov_b_1;
filter_retro.Pois{Ncom+index_extra}.t_bPois=-1; %We mark the time of birth with flag -1 (these will be integrated out)
filter_retro.Pois{Ncom+index_extra}.length_Pois=1;
filter_retro.Pois{Ncom+index_extra}.is_oos=1;




Ntracks=length(filter_upd.tracks);

if(Ntracks>0)
    
    for i=1:Ntracks
        
        Nhyp_i=length(filter_upd.tracks{i}.eB);
        
        %t_ini,t_b and length do not change
        t_bi=filter_upd.tracks{i}.t_b;
        
        
        if(t_bi>ko)
            %Trajectory is not present at OOS time. We do not do anything.
            
        elseif(t_bi==ko)
            %Scenario of non-present/present
            %We add an OOS state to all possible hypotheses (but we have
            %two hypotheses, one corresponding to non-extension and the
            %other corresponding to extension
            
            
            for j=1:Nhyp_i %We go through all hypotheses
                %We have two hypotheses, the first one is the same as
                %before (which is stored in variable meanB and covB. The new one
                %is in variable meanB2, covB2.
                mean_i=filter_retro.tracks{i}.meanB{j};
                cov_i=filter_retro.tracks{i}.covB{j};
                
                [mean_o,cov_o]=retro_np(mean_i,cov_i,Knp,Qnp,xb_prov_np,Nx);
                filter_retro.tracks{i}.meanB2{j}=mean_o;
                filter_retro.tracks{i}.covB2{j}=cov_o;
                
                %We need to do the past trajectories in the L-scan window
                Npast_means=length(filter_retro.tracks{i}.mean_past{j});
                
                for p=1:Npast_means
                    mean_i=filter_retro.tracks{i}.mean_past{j}{p};
                    cov_i=filter_retro.tracks{i}.cov_past{j}{p};
                    
                    [mean_o,cov_o]=retro_np(mean_i,cov_i,Knp,Qnp,xb_prov_np,Nx);
                    
                    filter_retro.tracks{i}.mean_past2{j}{p}=mean_o;
                    filter_retro.tracks{i}.cov_past2{j}{p}=cov_o;
                    
                end
                
            end
            
            
            
        else
            %The trajectory was born before ko. We have several cases for
            %each possible time of death
            
            filter_retro.tracks{i}.t_f=zeros(Nhyp_i,1);
            
            for j=1:Nhyp_i %We go through all hypotheses
                %We need to calculate final time step
                mean_i=filter_retro.tracks{i}.meanB{j};
                
                %Final time step
                t_f=t_bi+length(mean_i)/Nx-1;
                
                filter_retro.tracks{i}.t_f(j)=t_f;
                
                %If t_f<ko-1 trajectory and past values do not change. We
                %do not have to do anything (as values are already
                %included in the output)
                
                
                if(t_f==ko-1)
                    %We have a case of present/non-present
                    cov_i=filter_retro.tracks{i}.covB{j};
                    [mean_o,cov_o]=retro_pn(mean_i,cov_i,F1,Q1,Nx);
                    %We add the new trajectories
                    filter_retro.tracks{i}.meanB2{j}=mean_o;
                    filter_retro.tracks{i}.covB2{j}=cov_o;
                    
                    %No do not need to do anything with past trajectories as we
                    %know they do not exist
                    
                elseif(t_f>=ko)
                    %We have a case of present-present
                    %In present/present, we do not obtain a mixture in the
                    %retrodiction
                    
                    cov_i=filter_retro.tracks{i}.covB{j};
                    [mean_o,cov_o]=retro_pp(mean_i,cov_i, Kpp,Qpp,Fpp_prov,Nx,ko,t_f);
                    
                    
                    filter_retro.tracks{i}.meanB{j}=mean_o;
                    filter_retro.tracks{i}.covB{j}=cov_o;
                    
                    %We continue with past trajectories
                    Npast_means=length(filter_retro.tracks{i}.mean_past{j});
                    
                    for p=1:Npast_means
                        mean_i=filter_retro.tracks{i}.mean_past{j}{p};
                        t_f=t_bi+length(mean_i)/Nx-1;
                        
                        
                        filter_retro.tracks{i}.t_f_past{j}{p}=t_f;
                        
                        
                        if(p<Lscan-1)
                            
                            cov_i=filter_retro.tracks{i}.cov_past{j}{p};
  
                            if(t_f==ko-1)
                                %We have a case of present/non-present
                                [mean_o,cov_o]=retro_pn(mean_i,cov_i,F1,Q1,Nx);
                                
                                %We add the new trajectories
                                filter_retro.tracks{i}.mean_past2{j}{p}=mean_o;
                                filter_retro.tracks{i}.cov_past2{j}{p}=cov_o;
                                
                            elseif(t_f>=ko)
                                
                                
                                %We have a case of present-present
                                filter_retro.tracks{i}.is_oos_past{j}{p}=1;
                                [mean_o,cov_o]=retro_pp(mean_i,cov_i, Kpp,Qpp,Fpp_prov,Nx,ko,t_f);
                                
                                %We overwrite the trajectories with the state
                                %at OOS time
                                filter_retro.tracks{i}.mean_past{j}{p}=mean_o;
                                filter_retro.tracks{i}.cov_past{j}{p}=cov_o;
                                
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
            end
        end
        
    end
    
end
end

function  [mean_o,cov_o]=retro_np(mean_i,cov_i,Knp,Qnp,xb_prov_np,Nx)
%Function to extend with OOS state taking place just before
%(non-present/present hypothesis)

F_mode=zeros(Nx,size(cov_i,1));

F_mode(1:Nx,1:Nx)=Knp;
cross=F_mode*cov_i;
cov_o=[cov_i,cross';cross,F_mode*cov_i*F_mode'+Qnp];

mean_p=mean_i(end-size(cov_i,1)+1:end);
mean_o=[mean_i;...
    xb_prov_np+F_mode*mean_p];

end


function [mean_o,cov_o]=retro_pn(mean_i,cov_i,F1,Q1,Nx)
%Function to extend with OOS state
%(present/non-present hypothesis)
F_mode=zeros(Nx,size(cov_i,1));
F_mode(:,end-Nx+1:end)=F1;

cross=F_mode*cov_i;
cov_o=[cov_i,cross';cross,F_mode*cov_i*F_mode'+Q1];

mean_p=mean_i(end-Nx+1:end);
mean_o=[mean_i;...
    F1*mean_p];

end



function [mean_o,cov_o]=retro_pp(mean_i,cov_i, Kpp,Qpp,Fpp_prov,Nx,ko,t_f)
%Function to extend with OOS state
%(present/present hypothesis)

length=size(cov_i,1);
F_mode=zeros(Nx,length);

%These are the positions to fill out the transition matrix operating on the
%trajectory
index_top=length/Nx-(t_f-ko);
index_bottom=index_top-1;

F_mode(:,Nx*(index_top-1)+1:Nx*index_top)=Kpp;
F_mode(:,Nx*(index_bottom-1)+1:Nx*index_bottom)=Fpp_prov;

cross=F_mode*cov_i;
cov_o=[cov_i,cross';cross,F_mode*cov_i*F_mode'+Qpp];

mean_p=mean_i(end-size(cov_i,1)+1:end);

mean_o=[mean_i;...
    F_mode*mean_p];

end