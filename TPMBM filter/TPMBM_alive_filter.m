%Demo with the Gaussian implementation of the trajectory Poisson multi-Bernoulli
%mixture (TPMBM) filter for alive trajectories


%References
% K. Granström, L. Svensson, Y. Xia, J. Williams and Á. F. García-Femández, "Poisson Multi-Bernoulli Mixture Trackers: Continuity Through Random Finite Sets of Trajectories," 2018 21st International Conference on Information Fusion (FUSION), Cambridge, 2018
% This implementation was used in Á. F. García-Fernández, L. Svensson, J. L.
% Williams, Y. Xia, K. Granström, "Trajectory Poisson multi-Bernoulli
% filters",IEEE Transactions on Signal Processing, 2020.


%Copyright (c) 2020, Angel F. Garcia-Fernandez
%All rights reserved.


clear
addpath('..\GOSPA code')
addpath('..\Assignment')
addpath('..\PMBM filter')
addpath('..\Trajectory metric')
addpath('..\Trajectory errors')


rand('seed',9)
randn('seed',9)
ScenarioWilliams15;


plot_figures=0;


Nx=4; %Single target state dimension
%Filter
T_pruning=0.0001;%Threshold for pruning multi-Bernoulli mixtures weights
T_pruningPois=10^(-5); %Threshold for pruning PHD of the Poisson component
Nhyp_max=200;  %Maximum number of hypotheses (MBM components)
gating_threshold=20; %Threshold for gating
existence_threshold=0.00001; %Existence threshold: Bernoulli components with existence below this threshold are removed


type_estimator=1; %Choose Estimator 1, 2 or 3 as defined in the paper
existence_estimation_threshold1=0.4; %Only for esimator 1

%GOSPA errors for the estimator
squared_gospa_t_tot=zeros(1,Nsteps);
squared_gospa_loc_t_tot=zeros(1,Nsteps); %Localisation error
squared_gospa_false_t_tot=zeros(1,Nsteps); %False target error
squared_gospa_mis_t_tot=zeros(1,Nsteps); %Misdetection error


%LP_metric errors for the estimator
squared_LP_metric_t_tot=zeros(1,Nsteps);
squared_LP_metric_loc_t_tot=zeros(1,Nsteps); %Localisation error
squared_LP_metric_fal_t_tot=zeros(1,Nsteps); %False target error
squared_LP_metric_mis_t_tot=zeros(1,Nsteps); %Misdetection error
squared_LP_metric_switch_t_tot=zeros(1,Nsteps); %Track switching error



%Parameters for the trajectory version
Lscan=5;



rand('seed',9)
randn('seed',9)

%We go through all Monte Carlo runs

for i=1:Nmc
    tic
    
    %     filter_pred.weightPois=lambda0;
    %     filter_pred.meanPois=means_b;
    %     filter_pred.covPois=covs_b;
    
    filter_pred=cell(0,1);
    
    filter_pred.Pois{1}.weightPois=lambda0;
    filter_pred.Pois{1}.meanPois=means_b;
    filter_pred.Pois{1}.covPois=covs_b;
    filter_pred.Pois{1}.t_bPois=1; %t_birth
    filter_pred.Pois{1}.length_Pois=1; %length of trajectory
    
    
    
    
    filter_pred.tracks=cell(0,1);
    filter_pred.globHyp=[];
    filter_pred.globHypWeight=[];
    N_hypotheses_t=zeros(1,Nmc);
    
    
    
    
    %Simulate measurements
    for k=1:Nsteps
        z=CreateMeasurement(X_truth(:,k),t_birth,t_death,p_d,l_clutter,Area,k,H,chol_R,Nx);
        z_t{k}=z;
    end
    
    %Perform filtering
    
    for k=1:Nsteps
        
        %Update
        z=z_t{k};
        
        
        filter_upd=TPMBM_update(filter_pred,z,H,R,p_d,k,gating_threshold,intensity_clutter,Nhyp_max,Lscan);
        
        %State estimation
        switch type_estimator
            case 1
                [X_estimate,t_b_estimate,length_estimate]=TPMBM_estimate1(filter_upd,existence_estimation_threshold1);
            case 2
                [X_estimate,t_b_estimate,length_estimate]=TPMBM_estimate2(filter_upd);
            case 3
                [X_estimate,t_b_estimate,length_estimate]=TPMBM_estimate3(filter_upd);
        end
        
        
        
        
        %Computation of the squared LP metric error
        [squared_LP_metric, LP_metric_loc, LP_metric_miss, LP_metric_fal, LP_metric_switch]=ComputeLP_metric_error(X_estimate,t_b_estimate, length_estimate,X_truth,t_birth,t_death,c_gospa,gamma_track_metric,k,Nx);
        squared_LP_metric_t_tot(k)=squared_LP_metric_t_tot(k)+squared_LP_metric;
        squared_LP_metric_loc_t_tot(k)=squared_LP_metric_loc_t_tot(k)+LP_metric_loc;
        squared_LP_metric_fal_t_tot(k)=squared_LP_metric_fal_t_tot(k)+LP_metric_fal;
        squared_LP_metric_mis_t_tot(k)=squared_LP_metric_mis_t_tot(k)+LP_metric_miss;
        squared_LP_metric_switch_t_tot(k)=squared_LP_metric_switch_t_tot(k)+LP_metric_switch;
        
        
        %Computation of squared GOSPA position error and its decomposition
        %Obtain ground truth state
        [squared_gospa,gospa_loc,gospa_mis,gospa_fal]=ComputeGOSPAerror_trajectory(X_estimate,t_b_estimate, length_estimate,X_truth,t_birth,t_death,c_gospa,k,Nx);
        %We sum the squared errors
        squared_gospa_t_tot(k)=squared_gospa_t_tot(k)+squared_gospa;
        squared_gospa_loc_t_tot(k)=squared_gospa_loc_t_tot(k)+gospa_loc;
        squared_gospa_false_t_tot(k)=squared_gospa_false_t_tot(k)+gospa_fal;
        squared_gospa_mis_t_tot(k)=squared_gospa_mis_t_tot(k)+gospa_mis;
        
        
        
        
        
        
        %Draw filter output
        %         display(k)
        %          DrawTrajectoryFilterEstimates(X_truth,t_birth,t_death,X_estimate,[100,200],[100,200],z,k)
        
        
        
        %Hypothesis reduction, pruning,normalisation
        filter_upd_pruned=TPMBM_pruning(filter_upd, T_pruning,T_pruningPois,Nhyp_max,existence_threshold);
        filter_upd=filter_upd_pruned;
        
        
        
        
        N_hypotheses_t(k)=length(filter_upd.globHypWeight);
        
        
        %Prediction
        filter_pred=TPMBM_prediction(filter_upd,F,Q,p_s,weights_b,means_b,covs_b,Lscan,k);
        
    end
    
    t=toc;
    display(['Completed iteration number ', num2str(i),' time ', num2str(t), ' sec'])
    
end

%save(['TPMBM_alive_pd',int2str(100*p_d),'_R',int2str(R(1,1)),'_clut',int2str(l_clutter),'_Nhyp_max',int2str(Nhyp_max),'_Lscan',int2str(Lscan)])

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




