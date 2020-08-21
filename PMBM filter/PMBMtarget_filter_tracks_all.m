%Demo with the implementation of the Poisson multi-Bernoulli mixture (PMBM)
%filter described in
%Á. F. García-Fernández, J. L. Williams, K. Granström and L. Svensson, "Poisson Multi-Bernoulli Mixture Filter: Direct Derivation and Implementation," in IEEE Transactions on Aerospace and Electronic Systems, vol. 54, no. 4, pp. 1883-1901, Aug. 2018.

%This version estimates the set of all trajectories sequentially by linking target state estimates of the same Bernoulli, which
%is created by a unique measurement. This information is available in the filter metadata, see variables filter_upd.tracks{i}.t_ini and
%filter_upd.tracks{i}.aHis{1}(1)

%These two variables identify a Bernoulli component i, which can be made
%explicit in the posterior via auxiliary variables, see
%Á. F. García-Fernández, L. Svensson, J. L. Williams, Y. Xia, K. Granström,
%"Trajectory Poisson multi-Bernoulli filters," IEEE Transactions on Signal
%Processing, 2020.





clear
addpath('..\GOSPA code')
addpath('..\Assignment')
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

if(type_estimator>1)
    error('sequential track formation only implemented for estimator 1')
end




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
%"Trajectory Poisson multi-Bernoulli filters," IEEE Transactions on Signal Processing, 2020.
%2) (error online=0) We calculate the final trajectory metric error at the last time step of the simulation
error_online=1;



%Explanation of the structures filter_pred and filter_upd

%filter_pred and filter_upd are two structures with the following fields
%filter_upd.weightPois (weights of the Poisson point process (PPP)
%meanPois (means of the Gaussian components of the PPP)
%filter_upd.covPois (covs of the Gaussian components of the PPP)
%filter_upd.globHyp (matrix the number of rows is the number of global hypotheses and
%the columns the number of Bernoulli components. Each row indicates the
%index of the single target hypothesis for each Bernoulli component in the
%global hypotheses. It is zero if this global hypothesis does not have a
%particular Bernoulli component.
%filter_upd.globHypWeight: vector with the weights of the global hypotheses
%filter_upd.tracks cell array with each of the Bernoulli components
%filter_upd.tracks{i} contains the information of the ith Bernoulli component
%filter_upd.tracks{i}.meanB: means of the i-th Bernoulli component for each
%single target hypothesis
%filter_upd.tracks{i}.covB: covariance matrices of the i-th Bernoulli component for each
%single target hypothesis
%filter_upd.tracks{i}.eB: existence probabilities of the i-th Bernoulli component for each
%single target hypothesis
%filter_upd.tracks{i}.t_ini: time of birth of this Bernoulli component
%filter_upd.tracks{i}.aHis{j}: history of data associations (indices to
%measurements or 0 if undetected) for the j-th single target hypothesis of
%the i-th Bernoulli component

rand('seed',9)
randn('seed',9)


%We go through all Monte Carlo runs

for i=1:Nmc
    tic
    
    filter_pred.weightPois=lambda0;
    filter_pred.meanPois=means_b;
    filter_pred.covPois=covs_b;
    
    filter_pred.tracks=cell(0,1);
    filter_pred.globHyp=[];
    filter_pred.globHypWeight=[];
    N_hypotheses_t=zeros(1,Nmc);
    
    
    %Initial trajectory estimates
    X_estimate=cell(0,1);
    t_b_estimate=[];
    length_estimate=[];
    tag_estimate=zeros(2,0); %Includes t_ini and index of first measurement
    
    %Simulate measurements
    for k=1:Nsteps
        z=CreateMeasurement(X_truth(:,k),t_birth,t_death,p_d,l_clutter,Area,k,H,chol_R,Nx);
        z_t{k}=z;
    end
    
    
    
    %Perform filtering
    
    for k=1:Nsteps
        
        %Update
        z=z_t{k};
        
        
        filter_upd=PoissonMBMtarget_update(filter_pred,z,H,R,p_d,k,gating_threshold,intensity_clutter,Nhyp_max);
        
        %Set of targets and set of trajectories estimation
        [X_estimate,t_b_estimate,length_estimate,tag_estimate]=PoissonMBMtarget_estimate1_tracks(filter_upd,existence_estimation_threshold1,X_estimate,t_b_estimate,length_estimate,tag_estimate,k,Nx);
        
        
        if(error_online==1)
            %Computation of the squared LP metric error
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
        %DrawTrajectoryFilterEstimates(X_truth,t_birth,t_death,X_estimate,[100,200],[100,200],z,k)
        
        
        %Hypothesis reduction, pruning,normalisation
        filter_upd_pruned=PoissonMBMtarget_pruning(filter_upd, T_pruning,T_pruningPois,Nhyp_max,existence_threshold);
        filter_upd=filter_upd_pruned;
        
        N_hypotheses_t(k)=length(filter_upd.globHypWeight);
        
        %Prediction
        filter_pred=PoissonMBMtarget_pred(filter_upd,F,Q,p_s,weights_b,means_b,covs_b);
        
    end
    
    t=toc;
    display(['Completed iteration number ', num2str(i),' time ', num2str(t), ' sec'])
    
end

%save(['PMBM_tracks_all_pd',int2str(100*p_d),'_R',int2str(R(1,1)),'_clut',int2str(l_clutter),'_Nhyp_max',int2str(Nhyp_max)])




if(plot_figures)
    figure(1)
    plot(1:Nsteps,sqrt(squared_LP_metric_t_tot/Nmc),'Linewidth',1.3)
    xlabel('Time step')
    ylabel('RMS TM error')
    grid on
    
    figure(2)
    plot(1:Nsteps,sqrt(squared_LP_metric_loc_t_tot/Nmc),'Linewidth',1.3)
    xlabel('Time step')
    ylabel('RMS TM localisation error')
    grid on
    
    
    figure(3)
    plot(1:Nsteps,sqrt(squared_LP_metric_mis_t_tot/Nmc),'Linewidth',1.3)
    xlabel('Time step')
    ylabel('RMS TM error (missed targets)')
    grid on
    
    figure(4)
    plot(1:Nsteps,sqrt(squared_LP_metric_fal_t_tot/Nmc),'Linewidth',1.3)
    xlabel('Time step')
    ylabel('RMS TM error (false targets)')
    grid on
    
    figure(5)
    plot(1:Nsteps,sqrt(squared_LP_metric_switch_t_tot/Nmc),'Linewidth',1.3)
    xlabel('Time step')
    ylabel('RMS TM error (switching cost)')
    grid on
    
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