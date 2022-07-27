%Demo with the implementation of the trajectory Poisson multi-Bernoulli
%mixture (TPMBM) filter for the set of alive trajectories with VMF direction-of-arrival measurements and
%Gaussian range based on the iterated posterior linearisation filter (IPLF) based on
%conditional moments.

%This version does not keep the variances of old states of the
%trajectory (beyond the Lscan window)



% Author Angel F. Garcia-Fernandez

%Relevant papers


% Á. F. García-Fernández, J. Ralph, P. Horridge, S. Maskell, "Gaussian
% trajectory PMBM filter with nonlinear measurements based on posterior
% linearisation" in Proceedings of the 25th International Conference on
% Information Fusion, 2022

% A. F García-Fernández,J. Ralph, P. Horridge, S. Maskell, 
% "A Gaussian filtering method for multi-target tracking with
% nonlinear/non-Gaussian measurements"
% IEEE Transactions on Aerospace and Electronic Systems, 2021, 57, 3539-3548

% K. Granström, L. Svensson, Y. Xia, J. Williams and Á. F. García-Femández, 
% "Poisson Multi-Bernoulli Mixture Trackers: Continuity Through Random Finite Sets of Trajectories," 
% 2018 21st International Conference on Information Fusion (FUSION), 2018, pp. 1-5

% Á. F. García-Fernández, L. Svensson, M. R. Morelande and S. Särkkä, "Posterior Linearization Filter: 
% Principles and Implementation Using Sigma Points," in
% IEEE Transactions on Signal Processing, vol. 63, no. 20, pp. 5561-5573,
% Oct.15, 2015.

% F. Tronarp, Á. F. García-Fernández and S. Särkkä, "Iterative Filtering and Smoothing 
% in Nonlinear and Non-Gaussian Systems Using Conditional Moments," 
% in IEEE Signal Processing Letters, vol. 25, no. 3, pp. 408-412, March 2018.

%This code makes use of function circ_vmrnd.m by Philipp Berens
%see
%P. Berens, "CircStat: A MATLAB toolbox for circular statistics," Journal
%of Statistical Software, vol. 31, pp. 1–21, Sep. 2009.

% The author would like to thank DSTL Grant no. 1000143726 for financial support.

clear
addpath('..\GOSPA code')
addpath('..\Assignment')
addpath('..\PMBM filter')
addpath('..\Trajectory metric')
addpath('..\Trajectory errors')
addpath('..\TPMBM filter')

rand('seed',9)
randn('seed',9)
ScenarioWilliams15_range_bearing;

gamma_track_metric=1; %Gamma parameter (track switching penalty) for the metric on sets of trajectories

%plot_figures=0 does not show figures with errors per time step
%plot_figures=1 shows figures with errors computed by the linear programming metric for sets of trajectories (and its decomposition) per time step
plot_figures=1;


Nx=4; %Single target state dimension
%Filter
T_pruning=0.0001;%Threshold for pruning multi-Bernoulli mixtures weights
T_pruningPois=10^(-5); %Threshold for pruning PHD of the Poisson component
Nhyp_max=200;  %Maximum number of hypotheses (MBM components)
gating_threshold=50; %Threshold for gating
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

%Calculations/Parameters for the VMF updates

W0=1/3; %Weight of sigma-point at the origin of the unscented transform
Nz=3;
Wn=(1-W0)/(2*Nx);
weights_sp=[W0,Wn*ones(1,2*Nx)]; %Weights sigma-points
%Ap_kappa (correction factor of moments)
p=2; %p=2 as we are on the circle
Ap_kappa=besseli(p/2,kappa)./besseli(p/2-1,kappa);
index=isnan(Ap_kappa);
Ap_kappa(index)=1-(p-1)./(2*kappa(index)); %For large kappa we apply Filip´s expansion

log_I_k=log(besseli(0,kappa));   %Value of I0(k) (necessary for normalising constants)
if(isinf(log_I_k))
    log_I_k=kappa-1/2*log(2*pi*kappa)+log(1+1/(8*kappa));          %We apply the asymptotic formula (Mardia book Eq. (10.3.5))
end
log_denom_likelihood=log_I_k+log(2*pi*R_range)/2; %Denominator of the true l(z|x)
log_denom_Gauss=Nz*log(2*pi)/2; %Denominator for some approximated Gaussian likelihoods

Nit_iplf=5; %Maximum number of IPLF iterations
threshold_iplf=10^(-2); %threshold to control IPLF iterations, based on KLD

%Measurement directed linearisation flag: 0: standard IPLF initiation, 1: measurement-directed initiation.
mesurement_dir_flag=1; 

%Flag for perfoming likelihood correction: 0 we do the Gaussian
%approximation to the likelihood with constant pd (taken at the prior mean)
%as in the UKF. 1: we perform the likelihood correction, see ""A Gaussian filtering method for multi-target tracking with
% nonlinear/non-Gaussian measurements"
likelihood_corr_flag=1;
R_kf=zeros(3);
R_kf(3,3)=R_range;



rand('seed',9)
randn('seed',9)

%We go through all Monte Carlo runs

for i=1:Nmc
    tic
   
    
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
        z=CreateMeasurement_range_bearing(X_truth(:,k),t_birth,t_death,p_d,l_clutter,k,Nx,x_s,kappa,R_range,range_min,delta_range);
        z_t{k}=z;
    end
    
    %Perform filtering
    
    for k=1:Nsteps
        
        %Update
        z=z_t{k};      
        
        filter_upd=TPMBM_update_range_bearing(filter_pred,z,p_d,k,gating_threshold,intensity_clutter,Nhyp_max,Lscan,...
            kappa,Ap_kappa,x_s,weights_sp,R_kf,Nx,Nz,Nit_iplf,log_denom_likelihood,threshold_iplf,log_denom_Gauss,mesurement_dir_flag,likelihood_corr_flag);
        
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
        
        %display(k)
        %Draw filter output
%          z_pos=[z(3,:).*z(1,:)+x_s(1);z(3,:).*z(2,:)+x_s(2)];        
%          DrawTrajectoryFilterEstimates(X_truth,t_birth,t_death,X_estimate,[100,200],[100,200],z_pos,k)
        
        
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


 save(['TPMBM_alive_rb_kappa',int2str(kappa),'_R_range',...
     int2str(R_range),'_clut',int2str(l_clutter),'_m_dir_flag', int2str(mesurement_dir_flag),'_l_corr_flag',int2str(likelihood_corr_flag), '_Nit_iplf', int2str(Nit_iplf),'_thre_iplf',int2str(10000*threshold_iplf),'_Lscan',int2str(Lscan)]);


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


