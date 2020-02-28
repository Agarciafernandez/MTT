%Implemenation of the Continuous-discrete Gaussian mixture PHD filter
%(CD-PHD) filter described in

%A. F. García-Fernández, S. Maskell, "Continuous-discrete multiple target filtering: PMBM, PHD and CPHD ?lter implementations," IEEE Transactions on Signal Processing, vol. 68, pp. 1300-1314, 2020.
% Open access version in https://www.liverpool.ac.uk/electrical-engineering-and-electronics/staff/angel-garcia-fernandez/publications/



%Performance is measured with the GOSPA metric (alpha=2) and its
%decomposition into localisation errors for properly detected targets,
%and costs for missed targets and false targets (only possible for alpha=2).

%A. S. Rahmathullah, Á. F. García-Fernández and L. Svensson, "Generalized optimal sub-pattern assignment metric," 2017 20th International Conference on Information Fusion (Fusion), Xi'an, 2017, pp. 1-8.
% Short video on GOSPA
% https://www.youtube.com/watch?v=M79GTTytvCM


%Copyright (c) 2020, Angel F. Garcia-Fernandez
%All rights reserved.


clear
addpath('..\GOSPA code')
addpath('..\PMBM filter')
addpath('..\PHD filter')

rand('seed',9)
randn('seed',9)

Scenario_continuous;


%Required for GM-PHD filter code
Ncom_max=30;
T_pruning=10^(-5);
T_merging=0.1;
Ncom_b=1; %Ncom_b must be equal to one in this implementation

%GOSPA errors for the estimator with highest hypothesis
squared_gospa_t_tot=zeros(1,Nsteps);
squared_gospa_loc_t_tot=zeros(1,Nsteps); %Localisation error
squared_gospa_fal_t_tot=zeros(1,Nsteps); %False target error
squared_gospa_mis_t_tot=zeros(1,Nsteps); %Misdetection error

%We obtain the appearance birth paramaters partitioned for position and
%velocity (required for the best Gaussian match)
mean_a_p=mean_a(1:2:end);
mean_a_v=mean_a(2:2:end);
P_a_pp=P_a(1:2:end,1:2:end);
P_a_vv=P_a(2:2:end,2:2:end);
P_a_pv=P_a(1:2:end,2:2:end);

%Probability of survival at each time step
p_s_k_t=exp(-mu*delta_tk_t);

%plot_figures=0 does not show figures with errors per time step
%plot_figures=1 shows figures with GOSPA errors (and its decomposition) per time step
plot_figures=1;


rand('seed',9)
randn('seed',9)

%We go through all Monte Carlo runs

for i=1:Nmc
    tic
    
    %Simulate measurements
    for k=1:Nsteps
        z=CreateMeasurement(X_truth(:,k),t_birth,t_death,p_d,l_clutter,Area,k,H,chol_R,Nx);
        z_t{k}=z;
    end
    
    [mean_b_k,cov_b_k]=BirthParameters_cd(q,delta_tk_t(1),mean_a_p,mean_a_v,P_a_pp,P_a_vv,P_a_pv,mu);
    
    %Initialisation
    Ncom_k=1;
    weights_k=zeros(1,Ncom_max);
    means_k=zeros(Nx,Ncom_max);
    covs_k=zeros(Nx,Nx,Ncom_max);
    logical_active_k=false(1,Ncom_max);
    
    weights_k(1)=lambda_new_born_t(1);
    means_k(:,1)=mean_b_k;
    covs_k(:,:,1)=cov_b_k;
    logical_active_k(1)=1;
    
    
    
    
    %Perform filtering
    
    for k=1:Nsteps
        
        %Update
        z=z_t{k};
        
        [weights_u, means_u,covs_u,Ncom_u,logical_actives_u]=GMPHD_filter_update(weights_k, means_k,covs_k,Ncom_k,logical_active_k,z,H,R,p_d,l_clutter,Area);
        
        %Estimation
        X_estimate=GMPHD_estimation2(weights_u, means_u);
        
        
        %Computation of squared GOSPA position error and its decomposition
        %Obtain ground truth state
        [squared_gospa,gospa_loc,gospa_mis,gospa_fal]=ComputeGOSPAerror(X_estimate,X_truth,t_birth,t_death,c_gospa,k);
        
        %We sum the squared errors
        squared_gospa_t_tot(k)=squared_gospa_t_tot(k)+squared_gospa;
        squared_gospa_loc_t_tot(k)=squared_gospa_loc_t_tot(k)+gospa_loc;
        squared_gospa_fal_t_tot(k)=squared_gospa_fal_t_tot(k)+gospa_fal;
        squared_gospa_mis_t_tot(k)=squared_gospa_mis_t_tot(k)+gospa_mis;
        
        %Draw filter output
        %      DrawFilterEstimates(X_truth,t_birth,t_death,X_estimate,[-100,1000],[-100,1000],z,k)
        %         pause
        
        %Pruning
        [weights_u, means_u,covs_u,Ncom_u,logical_actives_u]=GMPHD_filter_pruning(weights_u, means_u,covs_u,logical_actives_u,Ncom_max,T_pruning,T_merging);
        
        
        %Prediction
        if(k<Nsteps)
            %We calculate the parameters of the transition density (bith
            %model, F, Q, and p_s
            [mean_b_k,cov_b_k]=BirthParameters_cd(q,delta_tk_t(k+1),mean_a_p,mean_a_v,P_a_pp,P_a_vv,P_a_pv,mu);
            p_s_k=p_s_k_t(k+1);
            F_k=kron(eye(2),[1 delta_tk_t(k+1);0 1]);
            Q_k=q*kron(eye(2),[delta_tk_t(k+1)^3/3 delta_tk_t(k+1)^2/2; delta_tk_t(k+1)^2/2 delta_tk_t(k+1)]);
            
            
            %Prediction (survival targets)
            [weights_k, means_k,covs_k,Ncom_k,logical_active_k]=GMPHD_filter_prediction(weights_u, means_u,covs_u,Ncom_u,logical_actives_u,F_k,Q_k,p_s_k);
            %Prediction (we add the new bon targets)
            
            %Ncom_b must be equal to one in this implementation
            index_birth=find(logical_active_k==0);
            Ncom_k=Ncom_k+Ncom_b;
            logical_active_k(index_birth(1:Ncom_b))=1;
            means_k(:,index_birth(1:Ncom_b))=mean_b_k;
            covs_k(:,:,index_birth(1:Ncom_b))=cov_b_k;
            weights_k(index_birth(1:Ncom_b))=lambda_new_born_t(k+1);
            
            
        end
    end
    
    t=toc;
    display(['Completed iteration number ', num2str(i),' time ', num2str(t), ' sec'])
    
end

%Root mean square GOSPA errors at each time step
rms_gospa_t=sqrt(squared_gospa_t_tot/Nmc);
rms_gospa_loc_t=sqrt(squared_gospa_loc_t_tot/Nmc);
rms_gospa_fal_t=sqrt(squared_gospa_fal_t_tot/Nmc);
rms_gospa_mis_t=sqrt(squared_gospa_mis_t_tot/Nmc);


%Root mean square GOSPA errors across all time steps
rms_gospa_tot=sqrt(sum(squared_gospa_t_tot)/(Nmc*Nsteps))
rms_gospa_loc_tot=sqrt(sum(squared_gospa_loc_t_tot)/(Nmc*Nsteps))
rms_gospa_fal_tot=sqrt(sum(squared_gospa_fal_t_tot)/(Nmc*Nsteps))
rms_gospa_mis_tot=sqrt(sum(squared_gospa_mis_t_tot)/(Nmc*Nsteps))

%save(['PHD_cd_pd',int2str(100*p_d),'_R',int2str(R(1,1)),'_clut',int2str(l_clutter),'_lambda',int2str(100*lambda),'_mu',int2str(100*mu)])

if(plot_figures==1)
    
    figure(1)
    plot(1:Nsteps,rms_gospa_t,'blue','Linewidth',1.3)
    grid on
    xlabel('Time step')
    ylabel('RMS GOSPA error')
    
    figure(2)
    plot(1:Nsteps,rms_gospa_loc_t,'blue','Linewidth',1.3)
    grid on
    xlabel('Time step')
    ylabel('RMS GOSPA localisation error')
    
    figure(3)
    plot(1:Nsteps,rms_gospa_fal_t,'blue','Linewidth',1.3)
    grid on
    xlabel('Time step')
    ylabel('RMS GOSPA false target error')
    
    figure(4)
    plot(1:Nsteps,rms_gospa_mis_t,'blue','Linewidth',1.3)
    grid on
    xlabel('Time step')
    ylabel('RMS GOSPA missed target error')
    
end
