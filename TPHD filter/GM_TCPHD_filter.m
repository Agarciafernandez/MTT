%Implementation of the Gaussian mixture Trajectory CPHD filter
%A. F. García-Fernández and L. Svensson, “Trajectory PHD and CPHD filters”,
%IEEE Transactions on Signal Processing, vol. 67, no. 22, pp. 5702-5714,
%Nov. 2019.

%This filter is based on random finite sets of trajectories, see
% A. F. García-Fernández, L. Svensson, and M. R. Morelande, “Multiple target tracking based on sets of trajectories,” IEEE Transactions on Aerospace and Electronic Systems. 
%[Online]. Available: https://arxiv.org/abs/1605.08163

%Evaluation using the metric for sets of trajectories, which allows to decompose the
%error into localisation errors, cost for missed targets, false targets,
%and track switches
% A. S. Rahmathullah, A. F. García-Fernández, and L. Svensson, “A metric on the space of ?nite sets of trajectories for evaluation of multi-target tracking algorithms,” 2016.
%[Online]. Available: http://arxiv.org/abs/1605.01177

%Also, evaluation using the GOSPA metric (alpha=2), which allows to decompose the
%error into localisation errors, cost for missed targets and false targets,
% A. S. Rahmathullah, A. F. García-Fernández, and L. Svensson, “Generalized optimal sub-pattern assignment metric,” in 20th International Conference on Information Fusion, 2017.

%Author: Angel F Garcia-Fernandez



clear
addpath('..\GOSPA code')
addpath('..\Trajectory metric')
addpath('..\Trajectory errors')
addpath('..\PMBM filter')


rand('seed',9)
randn('seed',9)
Scenario_GMTPHD_4targets;

%Filter parameters
Ncom_max=30;
T_pruning=10^(-4);
T_absorption=4;


%Parameters of metric evaluation (GOSPA and trajectory metric)
c_gospa=10;
gamma_track_metric=1;


%plot_figures=0 does not show figures with errors per time step
%plot_figures=1 shows figures with GOSPA errors (and its decomposition) per time step
%plot_figures=2 shows figures with trajectory metric errors (and its decomposition) per time step
plot_figures=0;


%Parameter L in the TPHD paper. It controls the length of the updated window
Lscan=5;


Ncardinality_max=10; %Maximum number of targets considered in cardinality distribution
lambda_birth=sum(weights_b);
cardinality_birth=(lambda_birth.^(0:Ncardinality_max)*exp(-lambda_birth)./factorial(0:Ncardinality_max))';


%GOSPA errors and decomposition (alpha=2)
squared_gospa_t_tot=zeros(1,Nsteps);
squared_gospa_loc_t_tot=zeros(1,Nsteps); %Localisation error
squared_gospa_fal_t_tot=zeros(1,Nsteps); %False target error
squared_gospa_mis_t_tot=zeros(1,Nsteps); %Misdetection error

%LP_metric errors for the estimator
squared_LP_metric_t_tot=zeros(1,Nsteps);
squared_LP_metric_loc_t_tot=zeros(1,Nsteps); %Localisation error
squared_LP_metric_fal_t_tot=zeros(1,Nsteps); %False target error
squared_LP_metric_mis_t_tot=zeros(1,Nsteps); %Misdetection error
squared_LP_metric_switch_t_tot=zeros(1,Nsteps); %Track switching error

rand('seed',9)
randn('seed',9)


for i=1:Nmc
    tic
    Ncom_k=0; %Number of PHD components at time step k
    weights_k=zeros(1,Ncom_max); %PHD weights at time step k
    means_k=zeros(Lscan*Nx,Ncom_max); %PHD mean at time step k (within the L-scan window) for each PHD component
    covs_k=zeros(Lscan*Nx,Lscan*Nx,Ncom_max); %Covariance matrix at time step k within the L-scan window for each PHD component
    logical_actives_k=false(1,Ncom_max); %Boolean that indicates if the i-th PHD component is empty or not
    t_ini_k=zeros(1,Ncom_max); %Initial time step of each PHD component
    length_k=zeros(1,Ncom_max); %Trajectory length of each PHD component 
    means_k_old=zeros(Nsteps*Nx,Ncom_max); %Contains the old means (before the window) I arrange them already according to its time step.
      
    
    for k=1:Nsteps
            
        %Birth
        
        index_birth=find(logical_actives_k==0);        
        Ncom_k=Ncom_k+Ncom_b;
        logical_actives_k(index_birth(1:Ncom_b))=1;
        means_k(1:Nx,index_birth(1:Ncom_b))=means_b;%The state is [xk-t;xk-t+1;...;xk]
        covs_k(1:Nx,1:Nx,index_birth(1:Ncom_b))=covs_b;
        t_ini_k(index_birth(1:Ncom_b))=k;
        length_k(index_birth(1:Ncom_b))=1;
        
        if(k>1)
            weights_k(index_birth(1:Ncom_b))=weights_b;
            cardinality_pred=CPHD_cardinality_pred(cardinality_k,p_s,cardinality_birth);
        else
            %The expected number of targets at time step zero is lambda0
            weights_k(index_birth(1:Ncom_b))=lambda0*weights_b/sum(weights_b);
            cardinality_pred=(lambda0.^(0:Ncardinality_max)*exp(-lambda0)./factorial(0:Ncardinality_max))';          
        end
        
        %Measurement generation
        z=CreateMeasurement(X_truth(:,k),t_birth,t_death,p_d,l_clutter,Area,k,H,chol_R,Nx);
        
        %Update   
        [weights_u, means_u,covs_u,t_ini_u,length_u,Ncom_u,logical_actives_u,cardinality_u,means_k_old_u]=GMTCPHD_Lscan_filter_update(weights_k, means_k,...
            covs_k,t_ini_k,length_k,Ncom_k,logical_actives_k,z,H,R,p_d,l_clutter,Area,cardinality_pred,Lscan,means_k_old);
               
        
        %Pruning and absorption
        [weights_o, means_o,covs_o,t_ini_o,length_o, Ncom_u,logical_actives_u,means_k_old_u]=GMTPHD_filter_pruning_absorption(weights_u,means_u,covs_u,t_ini_u,length_u,logical_actives_u,T_pruning,Ncom_max,...
            T_absorption,means_k_old_u,Lscan);
        weights_u=weights_o;
        means_u=means_o;
        covs_u=covs_o;
        t_ini_u=t_ini_o;
        length_u=length_o;
        
           
        %After pruning, the updated weights, are renormalised so that
        %their sum matches the cardinality mean
        weights_u=weights_u/sum(weights_u)*sum(cardinality_u.*(0:length(cardinality_u)-1)');
        
        
        %Estimation
        [X_estimate,t_b_estimate,length_estimate]=GMTCPHD_estimation(weights_u, means_u,t_ini_u,length_u,Lscan,means_k_old_u,cardinality_u);
        
        %Computation of the squared LP metric error
        
        [squared_LP_metric, LP_metric_loc, LP_metric_mis, LP_metric_fal, LP_metric_switch]=ComputeLP_metric_error(X_estimate,t_b_estimate, length_estimate,X_truth,t_birth,t_death,c_gospa,gamma_track_metric,k,Nx);
        squared_LP_metric_t_tot(k)=squared_LP_metric_t_tot(k)+squared_LP_metric;
        squared_LP_metric_loc_t_tot(k)=squared_LP_metric_loc_t_tot(k)+LP_metric_loc;
        squared_LP_metric_fal_t_tot(k)=squared_LP_metric_fal_t_tot(k)+LP_metric_fal;
        squared_LP_metric_mis_t_tot(k)=squared_LP_metric_mis_t_tot(k)+LP_metric_mis;
        squared_LP_metric_switch_t_tot(k)=squared_LP_metric_switch_t_tot(k)+LP_metric_switch;
        
        %Computation of squared GOSPA position error and its decomposition
        %Obtain ground truth state
        [squared_gospa,gospa_loc,gospa_mis,gospa_fal]=ComputeGOSPAerror_trajectory(X_estimate,t_b_estimate, length_estimate,X_truth,t_birth,t_death,c_gospa,k,Nx);
        
        %We sum the squared errors
        squared_gospa_t_tot(k)=squared_gospa_t_tot(k)+squared_gospa;
        squared_gospa_loc_t_tot(k)=squared_gospa_loc_t_tot(k)+gospa_loc;
        squared_gospa_fal_t_tot(k)=squared_gospa_fal_t_tot(k)+gospa_fal;
        squared_gospa_mis_t_tot(k)=squared_gospa_mis_t_tot(k)+gospa_mis;      
        
        %Draw filter output
        %DrawTrajectoryFilterEstimates(X_truth,t_birth,t_death,X_estimate,[-50,700],[-50,700],z,k)
        %pause(0.3)
        
        %Prediction (without new born trajectories and without cardinality)
        if(k<Nsteps)
            [weights_k, means_k,covs_k,t_ini_k,length_k,Ncom_k,logical_actives_k,means_k_old]=GMTPHD_Lscan_filter_prediction(weights_u, means_u,covs_u,...
                t_ini_u,length_u,Ncom_u,logical_actives_u,F,Q,p_s,Lscan,means_k_old_u,k);
            cardinality_k=cardinality_u;     
        end      
    end
       
    t=toc;
    display(['Completed iteration number ', num2str(i),' time ', num2str(t), ' sec'])
    
end


%save(['TCPHD_pd',int2str(100*p_d),'_R',int2str(R(1,1)),'_clut',int2str(l_clutter),'_Lscan',int2str(Lscan)])


display('LP trajectory metric errors')
%RMS trajectory metric costs across all time steps (for alive trajectories)
rms_LP_metric_tot=sqrt(sum(squared_LP_metric_t_tot)/(Nmc*Nsteps))
rms_LP_metric_loc_tot=sqrt(sum(squared_LP_metric_loc_t_tot)/(Nmc*Nsteps))
rms_LP_metric_mis_tot=sqrt(sum(squared_LP_metric_mis_t_tot)/(Nmc*Nsteps))
rms_LP_metric_false_tot=sqrt(sum(squared_LP_metric_fal_t_tot)/(Nmc*Nsteps))
rms_LP_metric_switch_tot=sqrt(sum(squared_LP_metric_switch_t_tot)/(Nmc*Nsteps))

display('GOSPA metric errors')
%RMS GOSPA costs across all time steps (for alive trajectories)
rms_gospa_tot=sqrt(sum(sum(squared_gospa_t_tot))/(Nmc*Nsteps))
rms_gospa_loc_tot=sqrt(sum(sum(squared_gospa_loc_t_tot))/(Nmc*Nsteps))
rms_gospa_mis_tot=sqrt(sum(sum(squared_gospa_mis_t_tot))/(Nmc*Nsteps))
rms_gospa_fal_tot=sqrt(sum(sum(squared_gospa_fal_t_tot))/(Nmc*Nsteps))


if(plot_figures==1)
    
    %Plot of GOSPA error and its decomposition
    figure(1)
    plot(1:Nsteps,sqrt(squared_gospa_t_tot/Nmc),'Linewidth',1.3)
    xlabel('Time step')
    ylabel('RMS GOSPA error')
    grid on
    
    figure(2)
    plot(1:Nsteps,sqrt(squared_gospa_loc_t_tot/Nmc),'Linewidth',1.3)
    xlabel('Time step')
    ylabel('RMS localisation cost (GOSPA)')
    grid on
        
    figure(3)
    plot(1:Nsteps,sqrt(squared_gospa_mis_t_tot/Nmc),'Linewidth',1.3)
    xlabel('Time step')
    ylabel('RMS missed target cost (GOSPA)')
    grid on
    
    figure(4)
    plot(1:Nsteps,sqrt(squared_gospa_fal_t_tot/Nmc),'Linewidth',1.3)
    xlabel('Time step')
    ylabel('RMS false target cost (GOSPA)')
    grid on
    
elseif(plot_figures==2)
    %Plot of LP trajectory metric error and its decomposition
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

