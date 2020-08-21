%Demo with the implementation of the trajectory multi-Bernoulli mixture (TMBM)
%filter with L-scan implementation and Murty algorithm

%Y. Xia, K. Granström, L. Svensson, A. F. García-Fernández, and J. L. Wlliams, “Multi-scan implementation of the trajectory Poisson multi-
%Bernoulli mixture filter,” Journal of Advances in Information Fusion. Special Issue on Multiple Hypothesis Tracking., vol. 14, no. 2, pp. 213–235, Dec. 2019.

%Implementation used in
%A. F. García-Fernández, L. Svensson, J. L. Williams, Y. Xia, K. Granström,  “Trajectory multi-Bernoulli filters for multi-target tracking based on sets of trajectories” in 23rd International Conference on Information Fusion, 2020.



%Filter for alive trajectories

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
ScenarioWilliams15_MB;


plot_figures=0;


Nx=4; %Single target state dimension
%Filter
T_pruning=0.0001;%Threshold for pruning multi-Bernoulli mixtures weights. This value does not affect the TPMB implementation as there is only one component. It is required for code compatibility
T_pruningPois=10^(-5); %Threshold for pruning PHD of the Poisson component
Nhyp_max=200;  %Maximum number of hypotheses (MBM components)
gating_threshold=20; %Threshold for gating
existence_threshold=0.00001; %Existence threshold: Bernoulli components with existence below this threshold are removed

%existence_threshold=0.0001;



type_estimator=1; %Choose Estimator 1, 2 or 3 as defined in the paper
existence_estimation_threshold1=0.4; %Only for esimator 1

%existence_estimation_threshold1=0.5;

T_alive=0.0001;



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

%We consider two options to evaluate GOSPA-trajectory metric error.
%1) (error_online=1) We calculate the error of the estimated set of all trajectories at each time step (normalised
%by the current time window) and add the error across time, see for example, (74) in
%Á. F. García-Fernández, L. Svensson, J. L. Williams, Y. Xia, K. Granström,
%"Trajectory Poisson multi-Bernoulli filters," accepted in IEEE-TSP.
%2) (error online=0) We calculate the final trajectory metric error at the last time step of the simulation
error_online=1;

%Parameters for the trajectory version
Lscan=5;


rand('seed',9)
randn('seed',9)

%We go through all Monte Carlo runs

for i=1:Nmc
    tic
    
    
    
    filter_pred=cell(0,1);
    filter_pred.tracks=cell(0,1);
    filter_pred.globHyp=ones(1,Ncom_b);
    filter_pred.globHypWeight=1;
    
    %Each Bernoulli component of the birth model starts a new track
    for j=1:Ncom_b
        filter_pred.tracks{j}.meanB{1}=means_b(:,j);
        filter_pred.tracks{j}.covB{1}=P_ini;
        filter_pred.tracks{j}.eB=e_ini(j);
        filter_pred.tracks{j}.t_ini=1;
        filter_pred.tracks{j}.aHis{1}=j;
        filter_pred.tracks{j}.t_b=1;
        filter_pred.tracks{j}.length=1;
        filter_pred.tracks{j}.prob_length{1}=1;
        filter_pred.tracks{j}.mean_past{1}=cell(0,1);
        
    end
    
    
    
    
    
    %Simulate measurements
    for k=1:Nsteps
        z=CreateMeasurement(X_truth(:,k),t_birth,t_death,p_d,l_clutter,Area,k,H,chol_R,Nx);
        z_t{k}=z;
    end
    
    %Perform filtering
    
    for k=1:Nsteps
        
        %Update
        z=z_t{k};
        
        
        
        filter_upd=TMBM_all_update(filter_pred,z,H,R,p_d,k,gating_threshold,intensity_clutter,Nhyp_max,Lscan,T_alive);
        
        
        %State estimation
        [X_estimate,t_b_estimate,length_estimate]=TPMBM_all_estimate1(filter_upd,existence_estimation_threshold1,Nx);
        
        
        if(error_online)
            
            %Error computation
            [squared_LP_metric, LP_metric_loc, LP_metric_miss, LP_metric_fal, LP_metric_switch]=ComputeLP_metric_all_error(X_estimate,t_b_estimate, length_estimate,X_truth,t_birth,t_death,c_gospa,gamma_track_metric,k,Nx);
            squared_LP_metric_t_tot(k)=squared_LP_metric_t_tot(k)+squared_LP_metric;
            squared_LP_metric_loc_t_tot(k)=squared_LP_metric_loc_t_tot(k)+LP_metric_loc;
            squared_LP_metric_fal_t_tot(k)=squared_LP_metric_fal_t_tot(k)+LP_metric_fal;
            squared_LP_metric_mis_t_tot(k)=squared_LP_metric_mis_t_tot(k)+LP_metric_miss;
            squared_LP_metric_switch_t_tot(k)=squared_LP_metric_switch_t_tot(k)+LP_metric_switch;
            %Computation of squared GOSPA position error and its decomposition
            %Obtain ground truth state
            [squared_gospa,gospa_loc,gospa_mis,gospa_fal]=ComputeGOSPAerror_all_trajectory(X_estimate,t_b_estimate, length_estimate,X_truth,t_birth,t_death,c_gospa,k,Nx);
            
            %We sum the squared errors
            squared_gospa_t_tot(k)=squared_gospa_t_tot(k)+squared_gospa;
            squared_gospa_loc_t_tot(k)=squared_gospa_loc_t_tot(k)+gospa_loc;
            squared_gospa_false_t_tot(k)=squared_gospa_false_t_tot(k)+gospa_fal;
            squared_gospa_mis_t_tot(k)=squared_gospa_mis_t_tot(k)+gospa_mis;
        else
            %We only calculate the errors at the last time step
            if(k==Nsteps)
                [squared_LP_metric, LP_metric_loc, LP_metric_miss, LP_metric_fal, LP_metric_switch]=ComputeLP_metric_all_error_final(X_estimate,t_b_estimate, length_estimate,X_truth,t_birth,t_death,c_gospa,gamma_track_metric,k,Nx);
                squared_LP_metric_t_tot=squared_LP_metric_t_tot+[zeros(Nsteps-1,1);squared_LP_metric];
                squared_LP_metric_loc_t_tot=squared_LP_metric_loc_t_tot+LP_metric_loc;
                squared_LP_metric_fal_t_tot=squared_LP_metric_fal_t_tot+LP_metric_fal;
                squared_LP_metric_mis_t_tot=squared_LP_metric_mis_t_tot+LP_metric_miss;
                squared_LP_metric_switch_t_tot=squared_LP_metric_switch_t_tot+[LP_metric_switch;0];
                
            end
            
        end
        
        %Draw filter output
        
        %DrawTrajectoryFilterEstimates(X_truth,t_birth,t_death,X_estimate,[50,250],[50,250],z,k)
        
        %Hypothesis reduction, pruning,normalisation
        filter_upd_pruned=TMBM_all_pruning(filter_upd, T_pruning,Nhyp_max,existence_threshold,T_alive);
        filter_upd=filter_upd_pruned;
        
        
        %Prediction
        
        filter_pred=TMBM_all_prediction(filter_upd,F,Q,p_s,means_b,P_ini,e_ini,Lscan,Ncom_b,k,T_alive,Nx);
        
        
    end
    
    t=toc;
    display(['Completed iteration number ', num2str(i),' time ', num2str(t), ' sec'])
    
end


%save(['TMBM_all_pd',int2str(100*p_d),'_R',int2str(R(1,1)),'_clut',int2str(l_clutter),'_Nhyp_max',int2str(Nhyp_max),'_Lscan',int2str(Lscan)])

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

if(error_online)
    
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
    
end
