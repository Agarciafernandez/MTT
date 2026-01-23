%Demo implementing several distributed Poisson multi-Bernoulli (PMB) filter variants using
%generalised covariance intersection (GCI), giving rise to these filters

% Distributed Poisson multi-Bernoulli mixture filter (DPMBM-GCI)
% Distributed PMB track oriented filter (DPMB-TO-GCI)
% Distributed PMB variational filter (DPMB-V-GCI)
% Distributed PMB global nearest neighbour filter (DPMB-GNN)

%The considered GCI fusion rule for PMBs and the characteristics of these filters are explained in the paper

% Á. F. García-Fernández and G. Battistelli, "Distributed Poisson Multi-Bernoulli Filtering via Generalized Covariance Intersection," in IEEE Transactions on Signal Processing, vol. 74, pp. 246-257, 2026, doi: 10.1109/TSP.2026.3651805.


%Copyright (c) 2025, Angel F. Garcia-Fernandez


clear
addpath('..\GOSPA code')
addpath('..\Assignment')
addpath('..\PMBM filter')
addpath('..\PHD filter') %This is required for the function GMPHD_filter_pruning.m, which is used to perform merging of PPP components


rand('seed',9)
randn('seed',9)

ScenarioWilliams15_distributed;
%Scenario_distributed;

Nagents=2; %This demo is only valid for 2 agents. It can be easily extended to multiple agents as explained in the paper.

if(Nagents>2)
    warning('This demo only works for two agents. See paper for the extension to multiple agents.')
end

Nsteps_fusion=10; %We do fusion every Nsteps_fusion


Nx=4; %Single target state dimension
%Filter
T_pruning=0.0001;%Threshold for pruning multi-Bernoulli mixtures weights
T_pruningPois=10^(-5); %Threshold for pruning PHD of the Poisson component
T_mergingPois=0.1; %Threshold for merging Poisson components
NcomPois_max=30; %Maximum number of PPP components

Nhyp_max=200;  %Maximum number of hypotheses (MBM components)
Nhyp_max_norm_prod=200; %Maximum number of hypotheses for PMBM of the normalised product of PMBs


gating_threshold=20; %Threshold for gating
existence_threshold=10^(-5); %Existence threshold: Bernoulli components with existence below this threshold are removed.


%existence_threshold=10^(-4) %For a higher number of sensors the threshold should be higher such that the algorithms is fast.


type_estimator=1; %Choose Estimator 1, 2 or 3 as defined in the paper
existence_estimation_threshold1=0.4; %Only for esimator 1

gating_threshold_fusion=20; %This is the gating threshold for fusion (normalised product)



%Flag to do the PMB projections
% (0) No PMB projection (we keep PMBM form)
% (1) TO-PMB projection
% (2) Variational PMB projection
pmb_proj_after_f=0; % Projection after the fusion step

%pmb_proj_after_f=2

pmb_proj_after_u=0; % Projection after the update step for each agent
% (0) No PMB projection (we keep PMBM form). If we use pmb_proj_after_u=0, we do a variational PMB projection before
%the fusion step, as this projection is what works best.
% (1) TO-PMB projection
% (2) Variational PMB projection



%Variational PMB (VPMB) projection parameters
max_iter_vpmb=10;
threshold_vpmb=0.1; %Threshold for convergence

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

rand('seed',9)
randn('seed',9)


omega=1/2;



%We go through all Monte Carlo runs

for i=1:Nmc
    tic

    %The information of each of the filters run by each agent is contained in this cell array
    filter_agents=cell(Nagents,1);
    for q=1:Nagents
        filter_agents{q}.filter_pred.weightPois=lambda0;
        filter_agents{q}.filter_pred.meanPois=means_b;
        filter_agents{q}.filter_pred.covPois=covs_b;

        filter_agents{q}.filter_pred.tracks=cell(0,1);
        filter_agents{q}.filter_pred.globHyp=[];
        filter_agents{q}.filter_pred.globHypWeight=[];

    end

    %Simulate measurements for each sensor
    measurement_agents=cell(Nagents,1);
    for q=1:Nagents
        for k=1:Nsteps
            z=CreateMeasurement(X_truth(:,k),t_birth,t_death,p_d,l_clutter,Area,k,H,chol_R,Nx);
            z_t{k}=z;
        end
        measurement_agents{q}=z_t;
    end


    %Perform filtering for each sensor

    for k=1:Nsteps

        %display(k)

        %Update for each sensor
        for q=1:Nagents
            %Set of measurements from sensor q at time step k
            z=measurement_agents{q}{k};

            filter_upd=PoissonMBMtarget_update(filter_agents{q}.filter_pred,z,H,R,p_d,k,gating_threshold,intensity_clutter,Nhyp_max);

            if(pmb_proj_after_u==1)

                %We project the PMBM updated density into a PMB density (option
                %1 track oriented form
                filter_upd_pmb=PMB_projection(filter_upd);
                filter_agents{q}.filter_upd=filter_upd_pmb;


            elseif(pmb_proj_after_u==2)
                %Option 2 variational form
                filter_upd_pruned=PMBMtarget_pruning_merging(filter_upd, T_pruning,T_pruningPois,Nhyp_max,existence_threshold,T_mergingPois,NcomPois_max);
                filter_upd_pmb=VPMB_projection(filter_upd_pruned,max_iter_vpmb,threshold_vpmb);

                %We do pruning to clean the PMB (some Bernoulli components may have
                %probability of existence=0 after the VPMB projection, and these must be removed before the next update)
                filter_upd_pmb=PMBMtarget_pruning_merging(filter_upd_pmb, T_pruning,T_pruningPois,Nhyp_max,existence_threshold,T_mergingPois,NcomPois_max);
                filter_agents{q}.filter_upd=filter_upd_pmb;

            else
                %Option 3: We keep PMBM form
                filter_agents{q}.filter_upd=filter_upd;
            end
        end



        if(rem(k,Nsteps_fusion)==0)
            %We perform fusion (done every Nsteps_fusion)

            if(pmb_proj_after_u==0)
                %If we are running PMBM filter, we do variational PMB
                %projection for each filter
                for q=1:Nagents
                    filter_upd_pruned=PMBMtarget_pruning_merging(filter_agents{q}.filter_upd, T_pruning,T_pruningPois,Nhyp_max,existence_threshold,T_mergingPois,NcomPois_max);
                    filter_upd_pmb=VPMB_projection(filter_upd_pruned,max_iter_vpmb,threshold_vpmb);

                    %We do pruning to clean the PMB (some Bernoulli components may have
                    %probability of existence=0 after the VPMB projection, and these must be removed before the next update)
                    filter_upd_pmb=PMBMtarget_pruning_merging(filter_upd_pmb, T_pruning,T_pruningPois,Nhyp_max,existence_threshold,T_mergingPois,NcomPois_max);
                    filter_agents{q}.filter_upd=filter_upd_pmb;
                end
            end

            filter_agents=PMB_GCI_Fusion_2_agents(filter_agents,omega,pmb_proj_after_f,gating_threshold_fusion,Nhyp_max_norm_prod,T_pruning,T_pruningPois,Nhyp_max,existence_threshold,T_mergingPois,NcomPois_max,max_iter_vpmb,threshold_vpmb);


        end

        %Plot multi-Bernoulli part of the PMB for both sensors (This function does
        %not work for PMBM
        % figure(6)
        % clf
        % hold on
        % DrawMB(filter_agents{1}.filter_upd,[0,300],[0,300],'r')
        % DrawMB(filter_agents{2}.filter_upd,[0,300],[0,300],'b')
        % hold off

        %State estimation
        X_estimate_agents=DPMBMStateEstimation(filter_agents,Nagents,type_estimator,existence_estimation_threshold1);

        %Computation of squared GOSPA position error and its
        %decomposition (the squared values are averaged over the number of
        %sensors at the end)

        for q=1:Nagents

            [squared_gospa,gospa_loc,gospa_mis,gospa_fal]=ComputeGOSPAerror(X_estimate_agents{q},X_truth,t_birth,t_death,c_gospa,k);

            %We sum the squared errors
            squared_gospa_t_tot(k)=squared_gospa_t_tot(k)+squared_gospa;
            squared_gospa_loc_t_tot(k)=squared_gospa_loc_t_tot(k)+gospa_loc;
            squared_gospa_false_t_tot(k)=squared_gospa_false_t_tot(k)+gospa_fal;
            squared_gospa_mis_t_tot(k)=squared_gospa_mis_t_tot(k)+gospa_mis;

        end


        %Draw distributed filters output
        % DrawDFilterEstimates(X_truth,t_birth,t_death,X_estimate_agents,[100,200],[100,200],measurement_agents,k,Nagents)
        % pause(0.5)


        %We perform hypothesis reduction, pruning, normalisation and prediction for each
        %filter
        for q=1:Nagents
            %Hypothesis reduction, pruning, merging of PPP components
            filter_upd=filter_agents{q}.filter_upd;
            filter_upd_pruned=PMBMtarget_pruning_merging(filter_upd, T_pruning,T_pruningPois,Nhyp_max,existence_threshold,T_mergingPois,NcomPois_max);
            filter_upd=filter_upd_pruned;
            filter_agents{q}.filter_upd=filter_upd;

            %Prediction
            filter_pred=PoissonMBMtarget_pred(filter_upd,F,Q,p_s,weights_b,means_b,covs_b);
            filter_agents{q}.filter_pred=filter_pred;

        end

    end
    % end
    t=toc;
    display(['Completed iteration number ', num2str(i),' time ', num2str(t), ' sec'])

end

%Root mean square GOSPA errors at each time step (also considering the
%number of agents)
rms_gospa_t=sqrt(squared_gospa_t_tot/(Nagents*Nmc));
rms_gospa_loc_t=sqrt(squared_gospa_loc_t_tot/(Nagents*Nmc));
rms_gospa_false_t=sqrt(squared_gospa_false_t_tot/(Nagents*Nmc));
rms_gospa_mis_t=sqrt(squared_gospa_mis_t_tot/(Nagents*Nmc));


%Root mean square GOSPA errors across all time steps (also considering the
%number of agents)
rms_gospa_tot=sqrt(sum(squared_gospa_t_tot)/(Nmc*Nagents*Nsteps))
rms_gospa_loc_tot=sqrt(sum(squared_gospa_loc_t_tot)/(Nmc*Nagents*Nsteps))
rms_gospa_false_tot=sqrt(sum(squared_gospa_false_t_tot)/(Nmc*Nagents*Nsteps))
rms_gospa_mis_tot=sqrt(sum(squared_gospa_mis_t_tot)/(Nmc*Nagents*Nsteps))


save(['DPMB_GCI_pd',int2str(100*p_d),'_R',int2str(R(1,1)),'_clut',int2str(l_clutter),'_Nhyp_max_norm_prod',int2str(Nhyp_max_norm_prod),'_pmb_proj_a_f',int2str(pmb_proj_after_f),'_pmb_proj_u',int2str(pmb_proj_after_u),'_Nsteps_fusion',int2str(Nsteps_fusion)])


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