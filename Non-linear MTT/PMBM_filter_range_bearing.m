%Demo with the implementation of the Poisson multi-Bernoulli mixture (PMBM)
%filter with VMF direction-of-arrival measurements and
%Gaussian range based on the iterated posterior linearisation filter (IPLF) based on
%conditional moments.

%Relevant references 

% A. F García-Fernández,J. Ralph, P. Horridge, S. Maskell,
% "A Gaussian filtering method for multi-target tracking with
% nonlinear/non-Gaussian measurements"
% IEEE Transactions on Aerospace and Electronic Systems, 2021, 57, 3539-3548

% J. L. Williams, "Marginal multi-Bernoulli filters: RFS derivation of MHT, JIPDA, and association-based member," 
% in IEEE Transactions on Aerospace and Electronic Systems, vol. 51, no. 3, pp. 1664-1687, July 2015,

% Á. F. García-Fernández, J. L. Williams, K. Granström and L. Svensson, "Poisson Multi-Bernoulli Mixture Filter:
% Direct Derivation and Implementation," in IEEE Transactions on Aerospace and Electronic Systems, vol. 54, no. 4, pp. 1883-1901, Aug. 2018.

% Á. F. García-Fernández, F. Tronarp and S. Särkkä, "Gaussian Target Tracking With Direction-of-Arrival von Mises–Fisher Measurements,"
% in IEEE Transactions on Signal Processing, vol. 67, no. 11, pp. 2960-2972, June, 2019.

% Á. F. García-Fernández, L. Svensson, M. R. Morelande and S. Särkkä, "Posterior Linearization Filter:
% Principles and Implementation Using Sigma Points," in
% IEEE Transactions on Signal Processing, vol. 63, no. 20, pp. 5561-5573,
% Oct.15, 2015.

% F. Tronarp, Á. F. García-Fernández and S. Särkkä, "Iterative Filtering and Smoothing 
% in Nonlinear and Non-Gaussian Systems Using Conditional Moments," 
% in IEEE Signal Processing Letters, vol. 25, no. 3, pp. 408-412, March 2018.

% The author would like to thank DSTL Grant no. 1000143726 for financial support.



%Author: Angel Garcia-Fernandez

clear
addpath('..\GOSPA code')
addpath('..\Assignment')
addpath('..\PMBM filter')

rand('seed',9)
randn('seed',9)
ScenarioWilliams15_range_bearing;


Nx=4; %Single target state dimension
%Filter
T_pruning=0.0001;%Threshold for pruning multi-Bernoulli mixtures weights
T_pruningPois=10^(-5); %Threshold for pruning PHD of the Poisson component
Nhyp_max=200;  %Maximum number of hypotheses (MBM components)
gating_threshold=50; %Threshold for gating

existence_threshold=10^(-4); %Existence threshold: Bernoulli components with existence below this threshold are removed

type_estimator=1; %Choose Estimator 1, 2 or 3 as defined in the paper
existence_estimation_threshold1=0.4; %Only for esimator 1

%GOSPA errors for the estimator with highest hypothesis
squared_gospa_t_tot=zeros(1,Nsteps);
squared_gospa_loc_t_tot=zeros(1,Nsteps); %Localisation error
squared_gospa_false_t_tot=zeros(1,Nsteps); %False target error
squared_gospa_mis_t_tot=zeros(1,Nsteps); %Misdetection error

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




Nit_iplf=5; %Maximum number of iterations
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

%plot_figures=0 does not show figures with errors per time step
%plot_figures=1 shows figures with GOSPA errors (and its decomposition) per time step
plot_figures=1;

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
    
    
    %Simulate measurements
    for k=1:Nsteps
        z=CreateMeasurement_range_bearing(X_truth(:,k),t_birth,t_death,p_d,l_clutter,k,Nx,x_s,kappa,R_range,range_min,delta_range);
        z_t{k}=z;
    end
    
    
    
    %Perform filtering
    
    for k=1:Nsteps
        
        %Update
        z=z_t{k};
        
        
        filter_upd=PoissonMBMtarget_update_range_bearing(filter_pred,z,p_d,k,gating_threshold,intensity_clutter,Nhyp_max,...
            kappa,Ap_kappa,x_s,weights_sp,R_kf,Nx,Nz,Nit_iplf,log_denom_likelihood,threshold_iplf,log_denom_Gauss,mesurement_dir_flag,likelihood_corr_flag);
        
        %State estimation
        switch type_estimator
            case 1
                X_estimate=PoissonMBMtarget_estimate1(filter_upd,existence_estimation_threshold1);
            case 2
                X_estimate=PoissonMBMtarget_estimate2(filter_upd);
            case 3
                X_estimate=PoissonMBMtarget_estimate3(filter_upd);
        end
        
        %Computation of squared GOSPA position error and its decomposition
        %Obtain ground truth state
        [squared_gospa,gospa_loc,gospa_mis,gospa_fal]=ComputeGOSPAerror(X_estimate,X_truth,t_birth,t_death,c_gospa,k);
        
     
        %We sum the squared errors
        squared_gospa_t_tot(k)=squared_gospa_t_tot(k)+squared_gospa;
        squared_gospa_loc_t_tot(k)=squared_gospa_loc_t_tot(k)+gospa_loc;
        squared_gospa_false_t_tot(k)=squared_gospa_false_t_tot(k)+gospa_fal;
        squared_gospa_mis_t_tot(k)=squared_gospa_mis_t_tot(k)+gospa_mis;
        
        
        %Uncomment to draw filter output at each time step
        
        %z_pos=[z(3,:).*z(1,:)+x_s(1);z(3,:).*z(2,:)+x_s(2)];
        %DrawFilterEstimates(X_truth,t_birth,t_death,X_estimate,[100,200],[100,200],z_pos,k)
        
        
        
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

%save(['PMBM_rb_kappa',int2str(kappa),'_R_range',int2str(R_range),'_clut',int2str(l_clutter),'_m_dir_flag', int2str(mesurement_dir_flag),'_l_corr_flag',int2str(likelihood_corr_flag), '_Nit_iplf', int2str(Nit_iplf),'_thre_iplf',int2str(10000*threshold_iplf)]);


%Root mean square GOSPA errors at each time step
rms_gospa_t=sqrt(squared_gospa_t_tot/Nmc);
rms_gospa_loc_t=sqrt(squared_gospa_loc_t_tot/Nmc);
rms_gospa_false_t=sqrt(squared_gospa_false_t_tot/Nmc);
rms_gospa_mis_t=sqrt(squared_gospa_mis_t_tot/Nmc);


%Root mean square GOSPA errors across all time steps
rms_gospa_tot=sqrt(sum(squared_gospa_t_tot)/(Nmc*Nsteps))
rms_gospa_loc_tot=sqrt(sum(squared_gospa_loc_t_tot)/(Nmc*Nsteps))
rms_gospa_false_tot=sqrt(sum(squared_gospa_false_t_tot)/(Nmc*Nsteps))
rms_gospa_mis_tot=sqrt(sum(squared_gospa_mis_t_tot)/(Nmc*Nsteps))


if(plot_figures)
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
    plot(1:Nsteps,rms_gospa_false_t,'blue','Linewidth',1.3)
    grid on
    xlabel('Time step')
    ylabel('RMS GOSPA false target error')
    
    figure(4)
    plot(1:Nsteps,rms_gospa_mis_t,'blue','Linewidth',1.3)
    grid on
    xlabel('Time step')
    ylabel('RMS GOSPA missed target error')
end
