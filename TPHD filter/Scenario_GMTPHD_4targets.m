%Scenario for the simulations

%Loads X_truth, t_birth, t_death, Nsteps
%X_truth,t_birth and t_death represents the set of ground truth
%trajectories, with their time of births and deaths
%Nsteps is the number of steps considered in the simulation
load('Sce_GMTPHD_4targets')

Nx=4; %State dimension

%Process model
T=0.5;
F=[1 T 0 0;0 1 0 0; 0 0 1 T; 0 0 0 1];
sigmaU=1.8;
Q=(sigmaU)^2*[T^3/3 T^2/2 0 0; T^2/2 T 0 0;0 0 T^3/3 T^2/2; 0 0 T^2/2 T];
p_s=0.99; %Probability of survival



%Measurement model (position measurements)
Area=[2000,2000]; %Area of interest
l_clutter=50; %Mean number of clutter measurements (Poisson cardinality and uniformly distributed in the surveillance area)
p_d=0.9; %Probability of detection
H=[1 0 0 0; 0 0 1 0];
R=diag([4 4]);
chol_R=chol(R)';

%Birth model
P_ini=diag([225 100 225 100]); %Covariance matrix of the birth intensity components
Ncom_b=3; %Number of components in the intensity of the birth mixture
weights_b=[0.1,0.1,0.1]; %Weights of the birth intensity mixture

%Expected number of targets at time step 0
lambda0=sum(weights_b);

%Means and covariance of the birth model
means_b=zeros(Nx,Ncom_b);
covs_b=zeros(Nx,Nx,Ncom_b);

covs_b(:,:,1)=P_ini;
covs_b(:,:,2)=P_ini;
covs_b(:,:,3)=P_ini;


means_b(:,1)=[85,0,140,0]'; 
means_b(:,2)=[-5,0,220,0]';
means_b(:,3)=[7,0,50,0]';





Nmc=100; %Number of Monte Carlo runs

