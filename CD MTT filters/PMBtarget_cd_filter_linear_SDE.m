%Implementation of the Continuous-discrete Poisson multi-Bernoulli 
%(CD-PMB) filter described in

%A. F. García-Fernández, S. Sarkka, "Gaussian multi-target filtering with
%target dynamics driven by a stochastic differential equation", in IEEE
%Transactions on Signal Processing, 2025.

%This code uses an Ornstein-Uhlenbeck (OU) process

%Author: Ángel F. García-Fernández


clear
addpath('..\GOSPA code')
addpath('..\Assignment')
addpath('..\PMBM filter')

rand('seed',9)
randn('seed',9)

Scenario_continuous_linear_SDE;

%Scenario_continuous_linear_SDE_v2;


%Filter
T_pruning=0.0001; %Threshold for pruning multi-Bernoulli mixtures weights
T_pruningPois=10^(-5); %Threshold for pruning PHD of the Poisson component
Nhyp_max=200;  %Maximum number of hypotheses (MBM components)
gating_threshold=20; %Threshold for gating
existence_threshold=10^(-5); %Existence threshold: Bernoulli components with existence below this threshold are removed

type_estimator=1; %Choose Estimator 1, 2 or 3 as defined in the paper
existence_estimation_threshold1=0.4; %Only for esimator 1

%GOSPA errors for the estimator with highest hypothesis
squared_gospa_t_tot=zeros(1,Nsteps);
squared_gospa_loc_t_tot=zeros(1,Nsteps); %Localisation error
squared_gospa_fal_t_tot=zeros(1,Nsteps); %False target error
squared_gospa_mis_t_tot=zeros(1,Nsteps); %Misdetection error



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
    
    %Initialisation

    [mean_b_k,cov_b_k]=BirthParameters_linear_SDE(mean_a,P_a,A,u,Q_c,mu, delta_tk_t(1));
    
    
    filter_pred.weightPois=lambda_new_born_t(1);
    filter_pred.meanPois=mean_b_k;
    filter_pred.covPois=cov_b_k;
    
    filter_pred.tracks=cell(0,1);
    filter_pred.globHyp=[];
    filter_pred.globHypWeight=[];
    N_hypotheses_t=zeros(1,Nmc);
    
    %Perform filtering
    
    for k=1:Nsteps
        

        %Update
        z=z_t{k};
        
        
        filter_upd=PoissonMBMtarget_update(filter_pred,z,H,R,p_d,k,gating_threshold,intensity_clutter,Nhyp_max);
        
        %We project the PMBM updated density into a PMB density
        filter_upd_pmb=PMB_projection(filter_upd);      
        filter_upd=filter_upd_pmb;

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
        squared_gospa_fal_t_tot(k)=squared_gospa_fal_t_tot(k)+gospa_fal;
        squared_gospa_mis_t_tot(k)=squared_gospa_mis_t_tot(k)+gospa_mis;
        
        %Draw filter output
        %        DrawFilterEstimates(X_truth,t_birth,t_death,X_estimate,[0,Area(1)],[0,Area(2)],z,k)
        %         pause
        
        %Hypothesis reduction, pruning,normalisation
        filter_upd_pruned=PoissonMBMtarget_pruning(filter_upd, T_pruning,T_pruningPois,Nhyp_max,existence_threshold);
        filter_upd=filter_upd_pruned;
        
        N_hypotheses_t(k)=length(filter_upd.globHypWeight);
        
        %Prediction
        if(k<Nsteps)
            %We calculate the parameters of the transition density (birth
            %model, F, Q, and p_s
            [mean_b_k,cov_b_k]=BirthParameters_linear_SDE(mean_a,P_a,A,u,Q_c,mu, delta_tk_t(k+1));
            p_s_k=p_s_k_t(k+1);




            [F_k,b_k,Q_k]=CalculateTransition_linear_SDE(A,u,delta_tk_t(k+1),Q_c);
  

            filter_pred=PoissonMBMtarget_pred_offset(filter_upd,F_k,Q_k,p_s_k,lambda_new_born_t(k+1),mean_b_k,cov_b_k,b_k);
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

save(['PMB_cd_linear_pd',int2str(100*p_d),'_R',int2str(R(1,1)),'_clut',int2str(l_clutter),'_Nhyp_max',int2str(Nhyp_max),'_lambda',int2str(100*lambda),'_mu',int2str(100*mu)])

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
