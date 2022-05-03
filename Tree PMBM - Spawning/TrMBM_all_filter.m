%Demo with the implementation of the Tree multi-Bernoulli
%mixture (TrMBM) with tree trajectories for multiple target target with
%spawning proposed in

%A. F. García-Fernández, L. Svensson, "Tracking multiple spawning targets
%using Poisson multi-Bernoulli mixtures on sets of tree trajectories", IEEE
%Transactions on Signal Processing, vol. 70, pp. 1987-1999, 2022.

%This version of the filter has a threshold T_alive. If the
%probability of the trajectory Bernoulli being alive at the current time step (given that it exists)
%is below this threshold, we consider the probability to be zero.
%Each single trajectory hypothesis of each Bernoulli component has another variable
%prob_length=[prob_length_k,prob_length_k-1,...,prob_lengthk-...] that tells us
%the probability of having a certain length in the last time steps. In the hypotheses
%that correspond to detection this variable is =1, as the trajectory is
%alive with probability one (in this Bernoulli component)

%The version for all trajectories also has another cell mean_past that stores the means considering the
%hypotheses that the trajectory died at a previous time step. For example, mean_past{1}{1} (stores mean at k-1 for hypothesis 1,) , mean_past{1}{2} (stores mean at k-1 for hypothesis 2).... that stores the
%previous means

%Author: A. F. García-Fernández


clear
addpath('..\GOSPA code')
addpath('..\Assignment')
addpath('..\PMBM filter')
addpath('..\TPMBM filter')
addpath('..\Trajectory metric')
addpath('..\Trajectory errors')


rand('seed',9)
randn('seed',9)
Scenario_spawning;

%Multi-Bernoulli birth parameters
Ncom_b=1;
e_ini=weights_b*ones(1,Ncom_b);
P_ini=covs_b;


plot_figures=0;


Nx=4; %Single target state dimension
%Filter
T_pruning=0.001;%Threshold for pruning multi-Bernoulli mixtures weights


T_pruningPois=10^(-4); %Threshold for pruning PHD of the Poisson component
Nhyp_max=100;  %Maximum number of hypotheses (MBM components)
gating_threshold=15; %Threshold for gating
existence_threshold=0.0001; %Existence threshold: Bernoulli components with existence below this threshold are removed

T_alive=0.0001;


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

%We consider two options to evaluate GOSPA-trajectory metric error.
%1) (error_online=1) We calculate the error of the estimated set of all trajectories at each time step (normalised
%by the current time window) and add the error across time, see for example, (74) in
%Á. F. García-Fernández, L. Svensson, J. L. Williams, Y. Xia, K. Granström,
%"Trajectory Poisson multi-Bernoulli filters," in IEEE-TSP.
%2) (error online=0) We calculate the final trajectory metric error at the last time step of the simulation

error_online=1;

%Parameters for the trajectory version
Lscan=1;

rand('seed',9)
randn('seed',9)

%We go through all Monte Carlo runs


for i=1:Nmc
    tic
    
 
    
    filter_pred=cell(0,1);
    filter_pred.tracks=cell(0,1);
    filter_pred.globHyp=ones(1,Ncom_b);
    filter_pred.globHypWeight=1;
    
    N_hypotheses_t=zeros(1,Nsteps);

    
    %Each Bernoulli component of the birth model starts a new track and
    %tree
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
        
        %Spawning variables
        filter_pred.tracks{j}.genealogy=1;
        
    end
    %Additional variables for spawning models
    %Tracks (which represent branches) have a variable called genealogy
    %(with the genealogy of the branch)
    %The variable tree index tell us the index of the tree that each branch
    %belongs to. The length of this variable is the number of
    %tracks/branches.
    filter_pred.tree_index=(1:Ncom_b)';
    filter_pred.N_tree=Ncom_b; %Number of Bernoulli trees, i.e., k*Ncomb_b for multi-Bernoulli birth model
    
    %Simulate measurements
    for k=1:Nsteps
        z=CreateMeasurement(X_truth(:,k),t_birth,t_death,p_d,l_clutter,Area,k,H,chol_R,Nx);
        z_t{k}=z;
    end
    
    %Perform filtering
    
    for k=1:Nsteps
        
        %Update
        z=z_t{k};
        

        
        filter_upd=TrMBM_all_update(filter_pred,z,H,R,p_d,k,gating_threshold,intensity_clutter,Nhyp_max,Lscan,T_alive);
        
        %State estimation
        switch type_estimator
            case 1
                [X_estimate,t_b_estimate,length_estimate]=TPMBM_all_estimate1(filter_upd,existence_estimation_threshold1,Nx);
            case 2
                %[X_estimate,t_b_estimate,length_estimate]=TPMBM_estimate2(filter_upd);
            case 3
                %[X_estimate,t_b_estimate,length_estimate]=TPMBM_estimate3(filter_upd);
        end
        
        
        if(error_online)
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
        
        
        
        
        %Draw filter output (two options, one for showing only branches,
        %the other for showing trees, uncomment required lines)
        
        %This function draws estimated branches (no genealogy information)
        %DrawTrajectoryFilterEstimates(X_truth,t_birth,t_death,X_estimate,[0,600],[0,400],z,k)
        
        
        %The next two lines draw estimated trees. First we need to estimate
        %the trees, which can be done by a modification of TPMBM_all_estimate1
        
%         [X_estimate,t_b_estimate,length_estimate,tree_index_estimate,genealogy_estimate]=TrPMBM_all_estimate1(filter_upd,existence_estimation_threshold1,Nx);
%         DrawTreeFilterEstimates(X_truth,t_birth,t_death,father_spawning,X_estimate,t_b_estimate,tree_index_estimate,genealogy_estimate,[0,600],[0,400],z,k)
        
        
        %Hypothesis reduction, pruning,normalisation
        filter_upd_pruned=TrMBM_all_pruning(filter_upd, T_pruning,Nhyp_max,existence_threshold,T_alive);
        filter_upd=filter_upd_pruned;
        
        N_hypotheses_t(k)=length(filter_upd.globHypWeight);
        
        
        %Prediction
        filter_pred=TrMBM_all_prediction(filter_upd,F,Q,p_s,means_b,P_ini,e_ini,Lscan,Ncom_b,k,T_alive,F1,F2,p_s1,p_s2,d_spawning);
        
    end
    
    t=toc;
    display(['Completed iteration number ', num2str(i),' time ', num2str(t), ' sec'])
    
end

save(['TrMBM_all_pd',int2str(100*p_d),'_R',int2str(R(1,1)),'_clut',int2str(l_clutter),'_Nhyp_max',int2str(Nhyp_max),'_Lscan',int2str(Lscan)])

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


