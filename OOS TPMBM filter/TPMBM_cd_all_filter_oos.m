%Continous-discrete trajectory Poisson multi-Bernoulli mixture (CD-TPMBM)
%filter with out-of-sequence (OOS) measurements described in

%A.F. Garcia-Fernandez, W. Yi "Continuous-discrete multiple target tracking with out-of-sequence measurements", 
% IEEE Transactions on Signal Processing, vol. 69, pp. 4699-4709, 2021.


%This filter keeps on information on the set of all (sampled) trajectories

%This version has a cell mean_past that stores the means considering the
%hypotheses that the trajectory died at a previous time step. For example, mean_past{1}{1} (stores mean at k-1 for hypothesis 1,) , mean_past{1}{2} (stores mean at k-1 for hypothesis 2).... that stores the
%previous means
%For the PPP, we only keep track of the alive trajectories

%In this version, we do not keep the variances of old states of the
%trajectory (beyond the Lscan window)

%We keep covariance matrices for trajectory hypotheses considering dead
%trajectories in the last L-scan window to perform OOS updates in this time
%window using cov_past in an analogous manner to mean_past.


%Author: Angel F. Garcia-Fernandez


clear
addpath('..\GOSPA code')
addpath('..\Assignment')
addpath('..\CD MTT filters')
addpath('..\PMBM filter')
addpath('..\Trajectory metric')
addpath('..\Trajectory errors')
addpath('..\TPMBM filter')


rand('seed',9)
randn('seed',9)
Scenario_oos;



plot_figures=0; %Set this to one to plot the figures
process_oos_flag=1; %If we set it to zero, we discard the oos measurements. It provides a baseline for comparison.


Nx=4; %Single target state dimension
%Filter
T_pruning=0.0001;%Threshold for pruning multi-Bernoulli mixtures weights
T_pruningPois=10^(-5); %Threshold for pruning PHD of the Poisson component

Nhyp_max=200;  %Maximum number of hypotheses (MBM components)
gating_threshold=20; %Threshold for gating
existence_threshold=0.00001; %Existence threshold: Bernoulli components with existence below this threshold are removed

T_alive=0.0001;



type_estimator=1; %Choose Estimator 1, 2 or 3 as defined in the paper
existence_estimation_threshold1=0.4; %Only for esimator 1

%GOSPA errors for the estimator
squared_gospa_t_tot=zeros(1,Nsteps);
squared_gospa_loc_t_tot=zeros(1,Nsteps); %Localisation error
squared_gospa_false_t_tot=zeros(1,Nsteps); %False target error
squared_gospa_mis_t_tot=zeros(1,Nsteps); %Misdetection error


%Linear programming (LP) trajectory metric errors for the estimator
squared_LP_metric_t_tot=zeros(1,Nsteps);
squared_LP_metric_loc_t_tot=zeros(1,Nsteps); %Localisation error
squared_LP_metric_fal_t_tot=zeros(1,Nsteps); %False target error
squared_LP_metric_mis_t_tot=zeros(1,Nsteps); %Misdetection error
squared_LP_metric_switch_t_tot=zeros(1,Nsteps); %Track switching error

%We obtain the appearance birth paramaters partitioned for position and
%velocity (required for the best Gaussian match). See
%Á. F. García-Fernández and S. Maskell, "Continuous-Discrete Multiple Target Filtering: PMBM, PHD and CPHD Filter Implementations," 
%in IEEE Transactions on Signal Processing, vol. 68, pp. 1300-1314, 2020.

mean_a_p=mean_a(1:2:end);
mean_a_v=mean_a(2:2:end);
P_a_pp=P_a(1:2:end,1:2:end);
P_a_vv=P_a(2:2:end,2:2:end);
P_a_pv=P_a(1:2:end,2:2:end);

%Probability of survival at each time step
p_s_k_t=exp(-mu*delta_tk_t);
p_s_k_t(~is_k_indices)=0; %For the oos-measurements, we do not use this p_s

%Intensity of new born targets at each time step (we overwrite the value
%with the oos values)
lambda_new_born_t=lambda/mu*(1-exp(-mu*delta_tk_t));
lambda_new_born_t(~is_k_indices)=0; %For the oos-measurements, we do not use this p_s




%L-scan Parameter for the trajectory implementation
Lscan=5;



rand('seed',9)
randn('seed',9)

%We go through all Monte Carlo runs

for i=1:Nmc
    tic
    
    %Simulate measurements
    for k=1:Nsteps
        k_oos=k_order(k); %We store them in the order we receive them
        z=CreateMeasurement(X_truth(:,k_oos),t_birth,t_death,p_d,l_clutter,Area,k,H,chol_R,Nx);
        z_t{k}=z;
    end
    [mean_b_k,cov_b_k]=BirthParameters_cd(q,delta_tk_t(1),mean_a_p,mean_a_v,P_a_pp,P_a_vv,P_a_pv,mu);
    
    %Initialisation
    
    filter_pred=cell(0,1);
    filter_pred.Pois{1}.weightPois=lambda_new_born_t(1);
    filter_pred.Pois{1}.meanPois=mean_b_k;
    filter_pred.Pois{1}.covPois=cov_b_k;
    filter_pred.Pois{1}.t_bPois=1; %t_birth
    filter_pred.Pois{1}.length_Pois=1; %length of trajectory
    
    filter_pred.tracks=cell(0,1);
    filter_pred.globHyp=[];
    filter_pred.globHypWeight=[];
    N_hypotheses_t=zeros(1,Nsteps);
    
    
    %Perform filtering
    
    for k=1:Nsteps
        
        
        %Update
        z=z_t{k};
        k_is=k_is_series(k); %Equivalent in-sequence time step.
        
        if(is_k_indices(k))
            %If Measurement is in sequence, we perform the standard update
            
            filter_upd=TPMBM_all_update_cov(filter_pred,z,H,R,p_d,k_is,gating_threshold,intensity_clutter,Nhyp_max,Lscan,T_alive);
            
            
        elseif(process_oos_flag)
            %process_oos_flag=1 indicates that we process the OOS
            %measurement
            
            to=tk(k);
            ko=find(tk_is>to,1);
            if(ko>=k_is-Lscan+2)
                %OOS is within Lscan window, we process the OOS measurement. Otherwise, we do not do anything
                
                %filter_upd represents the previous posterior
                filter_retro=TPMBM_all_retro(filter_upd,to,tk_is,k,k_is,q,mean_a_p,mean_a_v,P_a_pp,P_a_vv,P_a_pv,lambda,mu,Lscan);
                
                filter_retro_upd=TPMBM_all_retro_update(filter_retro,z,H,R,p_d,gating_threshold,intensity_clutter,Nhyp_max,Lscan);
                filter_upd=TPMBM_all_retro_marginal(filter_retro_upd,Nx,Lscan,k_is);
                
                
            end
            
        end
        
        %State estimation
        switch type_estimator
            case 1
                [X_estimate,t_b_estimate,length_estimate]=TPMBM_all_estimate1(filter_upd,existence_estimation_threshold1,Nx);
            case 2
                %[X_estimate,t_b_estimate,length_estimate]=TPMBM_estimate2(filter_upd);
            case 3
                %[X_estimate,t_b_estimate,length_estimate]=TPMBM_estimate3(filter_upd);
        end
        
        
        %Computation of the squared LP metric error
        [squared_LP_metric, LP_metric_loc, LP_metric_miss, LP_metric_fal, LP_metric_switch]=ComputeLP_metric_all_error(X_estimate,...
            t_b_estimate, length_estimate,X_truth_is,t_birth_is,t_death_is,c_gospa,gamma_track_metric,k_is,Nx);
        squared_LP_metric_t_tot(k)=squared_LP_metric_t_tot(k)+squared_LP_metric;
        squared_LP_metric_loc_t_tot(k)=squared_LP_metric_loc_t_tot(k)+LP_metric_loc;
        squared_LP_metric_fal_t_tot(k)=squared_LP_metric_fal_t_tot(k)+LP_metric_fal;
        squared_LP_metric_mis_t_tot(k)=squared_LP_metric_mis_t_tot(k)+LP_metric_miss;
        squared_LP_metric_switch_t_tot(k)=squared_LP_metric_switch_t_tot(k)+LP_metric_switch;
        
        %Computation of squared GOSPA position error and its decomposition
        %Obtain ground truth state
        [squared_gospa,gospa_loc,gospa_mis,gospa_fal]=ComputeGOSPAerror_all_trajectory(X_estimate,t_b_estimate,...
            length_estimate,X_truth_is,t_birth_is,t_death_is,c_gospa,k_is,Nx);
        
        %We sum the squared errors
        squared_gospa_t_tot(k)=squared_gospa_t_tot(k)+squared_gospa;
        squared_gospa_loc_t_tot(k)=squared_gospa_loc_t_tot(k)+gospa_loc;
        squared_gospa_false_t_tot(k)=squared_gospa_false_t_tot(k)+gospa_fal;
        squared_gospa_mis_t_tot(k)=squared_gospa_mis_t_tot(k)+gospa_mis;
        
        
        
        %Draw filter output
        
        %DrawTrajectoryFilterEstimates(X_truth_is,t_birth_is,t_death_is,X_estimate,[0,600],[0,400],z,k_is)
        
        
        %Hypothesis reduction, pruning,normalisation
        filter_upd_pruned=TPMBM_all_pruning_cov(filter_upd, T_pruning,T_pruningPois,Nhyp_max,existence_threshold,T_alive);
        filter_upd=filter_upd_pruned;
        
        
        
        
        N_hypotheses_t(k)=length(filter_upd.globHypWeight);

        
        %Prediction
        if(k<Nsteps)
            if(is_k_indices(k+1))
                %If the next measurement is in sequence, we perform prediction
                
                %We calculate the parameters of the transition density (bith
                %model, F, Q, and p_s
                [mean_b_k,cov_b_k]=BirthParameters_cd(q,delta_tk_t(k+1),mean_a_p,mean_a_v,P_a_pp,P_a_vv,P_a_pv,mu);
                p_s_k=p_s_k_t(k+1);
                F_k=kron(eye(2),[1 delta_tk_t(k+1);0 1]);
                Q_k=q*kron(eye(2),[delta_tk_t(k+1)^3/3 delta_tk_t(k+1)^2/2; delta_tk_t(k+1)^2/2 delta_tk_t(k+1)]);
                weights_b=lambda_new_born_t(k+1);
                %Below, we mark the new born trajectories with k_is(k)+1,as they are marked with the in-sequence measurements
                filter_pred=TPMBM_all_prediction_cov(filter_upd,F_k,Q_k,p_s_k,weights_b,mean_b_k,cov_b_k,Lscan,k_is,T_alive);
                
                
                
            end
            
        end
        
    end
    
    t=toc;
    display(['Completed iteration number ', num2str(i),' time ', num2str(t), ' sec'])
    
end



save(['TPMBM_oos_all_pd',int2str(100*p_d),'_R',int2str(R(1,1)),'_clut',int2str(l_clutter),'_Nhyp_max',int2str(Nhyp_max),'_Lscan',int2str(Lscan),'_oos',int2str(process_oos_flag),'_lambdaoos',int2str(lambda_oos)])

if(plot_figures)
    
    figure(1)
    plot(1:Nsteps,sqrt(squared_LP_metric_t_tot/Nmc),'blue','Linewidth',1.3)
    grid on
    xlabel('Time step')
    ylabel('RMS LP error')
    
    figure(2)
    plot(1:Nsteps,sqrt(squared_LP_metric_loc_t_tot/Nmc),'blue','Linewidth',1.3)
    grid on
    xlabel('Time step')
    ylabel('RMS LP error: localisation')
    
    
    figure(3)
    plot(1:Nsteps,sqrt(squared_LP_metric_fal_t_tot/Nmc),'blue','Linewidth',1.3)
    grid on
    xlabel('Time step')
    ylabel('RMS LP error: false targets')
    
    figure(5)
    plot(1:Nsteps,sqrt(squared_LP_metric_mis_t_tot/Nmc),'blue','Linewidth',1.3)
    grid on
    xlabel('Time step')
    ylabel('RMS LP error: missed targets')
    
    figure(6)
    plot(1:Nsteps,sqrt(squared_LP_metric_switch_t_tot/Nmc),'blue','Linewidth',1.3)
    grid on
    xlabel('Time step')
    ylabel('RMS LP error: track switches')
    
end

display('LP trajectory metric errors')
sqrt(sum(squared_LP_metric_t_tot)/(Nmc*Nsteps))
sqrt(sum(squared_LP_metric_loc_t_tot)/(Nmc*Nsteps))
sqrt(sum(squared_LP_metric_fal_t_tot)/(Nmc*Nsteps))
sqrt(sum(squared_LP_metric_mis_t_tot)/(Nmc*Nsteps))
sqrt(sum(squared_LP_metric_switch_t_tot)/(Nmc*Nsteps))

display('GOSPA metric errors')
sqrt(sum(squared_gospa_t_tot)/(Nmc*Nsteps))
sqrt(sum(squared_gospa_loc_t_tot)/(Nmc*Nsteps))
sqrt(sum(squared_gospa_false_t_tot)/(Nmc*Nsteps))
sqrt(sum(squared_gospa_mis_t_tot)/(Nmc*Nsteps))


