%Implemenation of the Continuous-discrete Gaussian mixture cardinality PHD filter
%(CD-CPHD) filter with dynamics following a non-linear SDE

%A. F. Garc�a-Fern�ndez, S. Sarkka, "Gaussian multi-target filtering with
%target dynamics driven by a stochastic differential equation", in IEEE
%Transactions on Signal Processing, 2025.

%For the specific example we use an nonlinear ODE with a re-entry tracking
%problem

%Author: �ngel F. Garc�a-Fern�ndez



clear
addpath('..\GOSPA code')
addpath('..\PMBM filter')
addpath('..\TPHD filter')
addpath('..\PHD filter')

rand('seed',9)
randn('seed',9)


Scenario_continuous_nonlinear_SDE;

%Scenario_continuous_nonlinear_SDE_v2;


%Required for GM-PHD filter code
Ncom_max=30;
T_pruning=10^(-5);
T_merging=0.1;
Ncom_b=1; %Ncom_b must be equal to one in this implementation
Ncardinality_max=15; %Maximum number of targets in cardinality

Ncardinality_max=40;

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

    [mean_b_k,cov_b_k]=BirthParameters_nonlinear_SDE_Taylor(mean_a,P_a, f_drift, F_drift,Q_c, mu, delta_tk_t(1));

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
    n_targets=(0:Ncardinality_max)';
    cardinality_pred=exp(-weights_k(1))*(weights_k(1)).^n_targets./factorial(n_targets); %Initial birth distribution is Poisson



    %Perform filtering

    for k=1:Nsteps

        %Update
        z=z_t{k};

        [weights_u, means_u,covs_u,Ncom_u,logical_active_u,cardinality_u]=GMCPHD_filter_update(weights_k, means_k,covs_k,Ncom_k,logical_active_k,z,H,R,p_d,l_clutter,Area,cardinality_pred);
        X_estimate=GMCPHD_estimation(weights_u, means_u,cardinality_u);


        %Computation of squared GOSPA position error and its decomposition
        %Obtain ground truth state
        [squared_gospa,gospa_loc,gospa_mis,gospa_fal]=ComputeGOSPAerror(X_estimate,X_truth,t_birth,t_death,c_gospa,k);

        %We sum the squared errors
        squared_gospa_t_tot(k)=squared_gospa_t_tot(k)+squared_gospa;
        squared_gospa_loc_t_tot(k)=squared_gospa_loc_t_tot(k)+gospa_loc;
        squared_gospa_fal_t_tot(k)=squared_gospa_fal_t_tot(k)+gospa_fal;
        squared_gospa_mis_t_tot(k)=squared_gospa_mis_t_tot(k)+gospa_mis;

        %Draw filter output
        %         DrawFilterEstimates(X_truth,t_birth,t_death,X_estimate,[-100,1000],[-100,1000],z,k)
        %         pause

        %Pruning
        [weights_u, means_u,covs_u,Ncom_u,logical_active_u]=GMPHD_filter_pruning(weights_u, means_u,covs_u,logical_active_u,Ncom_max,T_pruning,T_merging);

        %After pruning, the updated weights, are renormalised so that
        %their sum matches the cardinality mean
        weights_u=weights_u/sum(weights_u)*sum(cardinality_u.*(0:length(cardinality_u)-1)');

        %Prediction
        if(k<Nsteps)

            %We calculate the mean and covariance of the birth model
            [mean_b_k,cov_b_k]=BirthParameters_nonlinear_SDE_Taylor(mean_a,P_a, f_drift, F_drift,Q_c, mu, delta_tk_t(k+1));

            %Probability of survival
            p_s_k=p_s_k_t(k+1);

            %Prediction (survival targets)
            [weights_k, means_k,covs_k,Ncom_k,logical_active_k]=GMPHD_filter_prediction_nonlinear_SDE_Taylor(weights_u, means_u,covs_u,Ncom_u,logical_active_u,f_drift, F_drift,Q_c,p_s_k,delta_tk_t(k+1),Nx);


            %Prediction (we add the new bon targets) PHD
            %Ncom_b must be equal to one in this implementation
            index_birth=find(logical_active_k==0);
            Ncom_k=Ncom_k+Ncom_b;
            logical_active_k(index_birth(1:Ncom_b))=1;
            means_k(:,index_birth(1:Ncom_b))=mean_b_k;
            covs_k(:,:,index_birth(1:Ncom_b))=cov_b_k;
            weights_k(index_birth(1:Ncom_b))=lambda_new_born_t(k+1);

            %Prediction: cardinality
            cardinality_birth=(lambda_new_born_t(k+1).^(0:Ncardinality_max)*exp(-lambda_new_born_t(k+1))./factorial(0:Ncardinality_max))';
            cardinality_pred=CPHD_cardinality_pred(cardinality_u,p_s_k,cardinality_birth);



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

save(['CPHD_cd_nonlinear_Taylor_pd',int2str(100*p_d),'_R',int2str(R(1,1)),'_clut',int2str(l_clutter),'_lambda',int2str(100*lambda),'_mu',int2str(100*mu)])


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
