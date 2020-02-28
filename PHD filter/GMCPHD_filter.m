%Implementation of the Gaussian mixture CPHD filter in 
%B.-T. Vo, B.-N. Vo, and A. Cantoni, “Analytic implementations of the
%cardinalized probability hypothesis density filter,” IEEE Trans. Signal Process., vol. 55, no. 7, pp. 3553–3567, Jul. 2007.

%Performance is measured with the GOSPA metric (alpha=2) and its
%decomposition into localisation errors for properly detected targets,
%and costs for missed targets and false targets (only possible for alpha=2).
%A. S. Rahmathullah, Á. F. García-Fernández and L. Svensson, "Generalized optimal sub-pattern assignment metric," 2017 20th International Conference on Information Fusion (Fusion), Xi'an, 2017, pp. 1-8.

%Author: Angel Garcia-Fernandez

clear
addpath('..\GOSPA code')
addpath('..\PMBM filter')
addpath('..\TPHD filter')

rand('seed',9)
randn('seed',9)
ScenarioWilliams15;


%Required for GM-CPHD filter code

Ncom_max=30;
Nx=size(F,1);
T_pruning=10^(-5);
T_merging=0.1;
Ncardinality_max=10; %Maximum number of targets in cardinality


lambda_birth=sum(weights_b);
cardinality_birth=(lambda_birth.^(0:Ncardinality_max)*exp(-lambda_birth)./factorial(0:Ncardinality_max))';

%GOSPA errors for the estimator with highest hypothesis
squared_gospa_t_tot=zeros(1,Nsteps);
squared_gospa_loc_t_tot=zeros(1,Nsteps); %Localisation error
squared_gospa_false_t_tot=zeros(1,Nsteps); %False target error
squared_gospa_mis_t_tot=zeros(1,Nsteps); %Misdetection error


for i=1:Nmc
    tic
    Ncom_k=0;
    weights_k=zeros(1,Ncom_max);
    means_k=zeros(Nx,Ncom_max);
    covs_k=zeros(Nx,Nx,Ncom_max);
    logical_actives_k=false(1,Ncom_max);
    cardinality_k=[1;zeros(Ncardinality_max,1)];
        
     %Simulate measurements
    for k=1:Nsteps
        z=CreateMeasurement(X_truth(:,k),t_birth,t_death,p_d,l_clutter,Area,k,H,chol_R,Nx);
        z_t{k}=z;
    end
    
    for k=1:Nsteps
        
    
        %Birth (Poisson birth model)
        
        index_birth=find(logical_actives_k==0);
        
        
        Ncom_k=Ncom_k+Ncom_b;
        logical_actives_k(index_birth(1:Ncom_b))=1;
        means_k(:,index_birth(1:Ncom_b))=means_b;
        covs_k(:,:,index_birth(1:Ncom_b))=covs_b;
              
        if(k>1)
            weights_k(index_birth(1:Ncom_b))=weights_b;
            cardinality_pred=CPHD_cardinality_pred(cardinality_k,p_s,cardinality_birth);
        else
             %The expected number of targets at time step zero is lambda0
            weights_k(index_birth(1:Ncom_b))=lambda0*weights_b/sum(weights_b);
            cardinality_pred=(lambda0.^(0:Ncardinality_max)*exp(-lambda0)./factorial(0:Ncardinality_max))';  
        end
        
        
        
         %Measurement generation
        z=z_t{k};
        
        %Update
        [weights_u, means_u,covs_u,Ncom_u,logical_actives_u,cardinality_u]=GMCPHD_filter_update(weights_k, means_k,covs_k,Ncom_k,logical_actives_k,z,H,R,p_d,l_clutter,Area,cardinality_pred);
        
        
        %Estimation
        X_estimate=GMCPHD_estimation(weights_u, means_u,cardinality_u);
        
        %Computation of squared GOSPA position error and its decomposition
        %Obtain ground truth state
        [squared_gospa,gospa_loc,gospa_mis,gospa_fal]=ComputeGOSPAerror(X_estimate,X_truth,t_birth,t_death,c_gospa,k);
        
        %We sum the squared errors
        squared_gospa_t_tot(k)=squared_gospa_t_tot(k)+squared_gospa;
        squared_gospa_loc_t_tot(k)=squared_gospa_loc_t_tot(k)+gospa_loc;
        squared_gospa_false_t_tot(k)=squared_gospa_false_t_tot(k)+gospa_fal;
        squared_gospa_mis_t_tot(k)=squared_gospa_mis_t_tot(k)+gospa_mis;
        
        
        %Draw filter output
        %DrawFilterEstimates(X_truth,t_birth,t_death,X_estimate,[100,200],[100,200],z,k)
        %pause
         
        
        [weights_u, means_u,covs_u,Ncom_u,logical_actives_u]=GMPHD_filter_pruning(weights_u, means_u,covs_u,logical_actives_u,Ncom_max,T_pruning,T_merging);
        
        %Prediction (survival targets)
        
        [weights_k, means_k,covs_k,Ncom_k,logical_actives_k]=GMPHD_filter_prediction(weights_u, means_u,covs_u,Ncom_u,logical_actives_u,F,Q,p_s);
        cardinality_k=cardinality_u;
        
        
    end
    
  
    t=toc;
    display(['Completed iteration number ', num2str(i),' time ', num2str(t), ' sec'])


    
end

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
