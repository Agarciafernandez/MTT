% Demo with the implementation of the point-target Poisson multi-Bernoulli mixture (PMBM)
% filter with negative binomial clutter described in
% A. F. Garcia-Fernandez, Y. Xia, L. Svensson, "Poisson multi-Bernoulli
% mixture filter with general target-generated measurements and arbitrary
% clutter", IEEE Transactions on Signal Processing, vol. 71, 2023.


%Author: Angel F. Garcia-Fernandez

clear
addpath('..\GOSPA code')
addpath('..\Assignment')
addpath('..\PMBM filter')


rand('seed',9)
randn('seed',9)

Scenario_nb_clutter;

plot_figures=1; %To plot figures
nb_indicator=1; %nb_indicator=1 -> we use the PMB for negative binomial clutter. nb_indicator=0 -> we use the standard PMB with Gibbs sampling for data associations.

if(nb_indicator)
    %We precompute the log pdf of the clutter for different number of
    %clutter measurements
    card_max=400; % Maximum number of possible clutter measurements (If it is set too low, it may return an error)
    log_pdf_clutter=zeros(card_max+1,1);
    for m=1:length(log_pdf_clutter)
        log_pdf_clutter(m)=Calculate_log_clutter_pdf_nb(m-1,r_nb,p_nb,Area(1)*Area(2));
    end
else
    intensity_clutter=l_clutter/(Area(1)*Area(2)); %This should only be used in the standard PMBM implementation, not for the nb-PMBM filter
end


Nx=4; %Single target state dimension
%Filter
T_pruning=0.0001;%Threshold for pruning multi-Bernoulli mixtures weights
T_pruningPois=10^(-5); %Threshold for pruning PHD of the Poisson component
Nhyp_max=5000;  %Maximum number of global hypotheses (MBM components)


gating_threshold=20; %Threshold for gating
existence_threshold=0.00001; %Existence threshold: Bernoulli components with existence below this threshold are removed

type_estimator=3; %Choose Estimator 1, 2 or 3 as defined in the paper
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



%We go through all Monte Carlo runs

for i=1:Nmc

    tic

    filter_pred.weightPois=lambda0;
    filter_pred.meanPois=means_b;
    filter_pred.covPois=covs_b;

    filter_pred.tracks=cell(0,1);
    filter_pred.globHyp=[];
    filter_pred.globHypWeight=[];
    N_hypotheses_t=zeros(1,Nsteps);


    rand('seed',i)
    randn('seed',i)


    %Simulate measurements
    for k=1:Nsteps
        z=CreateMeasurement_nb_clutter(X_truth(:,k),t_birth,t_death,p_d,Area,k,H,chol_R,Nx,r_nb,p_nb);
        z_t{k}=z;
    end

    %Perform filtering

    for k=1:Nsteps


        %Update
        z=z_t{k};

        if(nb_indicator)
            filter_upd=PoissonMBMtarget_update_nb_clutter(filter_pred,z,H,R,p_d,k,gating_threshold,Nhyp_max,log_pdf_clutter);
        else
            filter_upd=PoissonMBMtarget_update_Gibbs(filter_pred,z,H,R,p_d,k,gating_threshold,intensity_clutter,Nhyp_max);
        end


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

        %Draw filter output
       % DrawFilterEstimates(X_truth,t_birth,t_death,X_estimate,[0,Area(1)],[0,Area(2)],z,k)
        %pause(0.5)



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

save(['PMBM_nb_pd',int2str(100*p_d),'_R',int2str(R(1,1)),'_clut',int2str(l_clutter), '_anb',int2str(a_nb),'_Nhyp_max',int2str(Nhyp_max),'_nb',int2str(nb_indicator),'_est',int2str(type_estimator)])


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
