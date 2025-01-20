rand('seed',9)
randn('seed',9)

%Simulation parameters
Nmc=100; %Number of Monte Carlo runs

%Measurement parameters

H=kron(eye(2),[1 0]);
R=4*eye(2);

R=2.5*10^(-3)*eye(2);

chol_R=chol(R)';
p_d=0.9;
Area=[900 6600];   
Nsteps=100; %Considered number of time steps in the simulation
l_clutter=10; %Average number of clutter measurements per time step



intensity_clutter=l_clutter/(Area(1)*Area(2)); %Intensity of the clutter

%Parameters of the continous time birth/death process

lambda=0.3;
mu=0.01;


%Mean and covariance matrix at time of appearance
% mean_a=[6500;-0.1;350;-0.1];
% P_a=diag([50 1 50 3].^2);

mean_a=[450;-0.1;6500;-1];
P_a=diag([100 3 50 1].^2);

%We generate the time intervals between different measurements (we take Nsteps measurements) 
%delta_tk (which is known but samples from an exponential distribution with parameter mu_delta)
mu_measurement=1; %Note that this is related to time through 1/mu_measurement
delta_tk_t=exprnd(1/mu_measurement,Nsteps,1);
tk=cumsum(delta_tk_t); %Time at which is measurement is obtained

%Uncomment to plot the time differences between measurements
% figure(7)
% plot(delta_tk_t)
% grid on
% xlabel('Time step')
% ylabel('Time interval \Deltat_{k} (s)')

%SDE parameters (We use the re-entry tracking problem)

Gm0=3.9860*10^(5);
beta0=0.6;
Psi_aero=0.7;
r0=6374;
h0=13.4;

%We set the drift function and its Jacobian
[f_drift,F_drift]=Set_drift_function(Gm0,beta0,Psi_aero,r0,h0);

step_size=0.01; %Step size of the Euler-Maruyama method for simulation

L=[0,0;...
    1,0;...
    0,0;...
    0, 1];


%Process noise intensity
q=0.0001;
Q_c=L*(q*eye(2))*L';


%Sampling of the multi-target system model

%We sample the number of new born targets at each time step
lambda_new_born_t=lambda/mu*(1-exp(-mu*delta_tk_t));
N_new_born_t=poissrnd(lambda_new_born_t);
N_targets_tot=sum(N_new_born_t);

%Life span for each target is exponential distrubuted
life_span_continuous=exprnd(1/mu,N_targets_tot,1); 

%We go through each target and generate its states, and vectors t_birth and t_death
Nx=4;

%X_truth matrix with the ground truth set of trajectories/targets
X_truth=zeros(Nx*N_targets_tot,Nsteps);
t_birth=zeros(1,N_targets_tot); %t_birth is the discrete time step at which the target is born
t_death=zeros(1,N_targets_tot); %t_death is the discrete time step the target is already dead

N_new_born_t_prov=N_new_born_t; 
time_aux=1; %Auxiliary variable
N_targets_alive_t=zeros(1,Nsteps);

%Exponential distribution
expon_dist=makedist('Exponential','mu',1/mu); %Note: makedist considers mu as the mean, so we need 1/mu

for i=1:N_targets_tot
    while(N_new_born_t_prov(time_aux)==0)
        time_aux=time_aux+1;
    end
    N_new_born_t_prov(time_aux)=N_new_born_t_prov(time_aux)-1;
    t_birth(i)=time_aux;
    t_death_continuous=life_span_continuous(i)+tk(time_aux);
    t_death_discrete=find(tk<t_death_continuous);
    t_death_discrete=t_death_discrete(end)+1; %We sum one
    t_death(i)=t_death_discrete;
    
    N_targets_alive_t(t_birth(i):t_death(i)-1)=N_targets_alive_t(t_birth(i):t_death(i)-1)+1;
    
    %Generation of the samples
    %Time difference between measurements at the current time step
    delta_tk_i=delta_tk_t(time_aux);
    
    %We truncate the exponential distribution according to delta_tk_i and
    %sample (this is the distribution of the time lags, see paper).
    expon_dist_trunc=truncate(expon_dist,0,delta_tk_i);    
    time_lag_i=random(expon_dist_trunc,1);
    

    %Sample of target state at appearing time
    X_a=mean_a+chol(P_a)'*randn(Nx,1);
    X_b=Euler_simulation(X_a,f_drift,L,q,time_lag_i,step_size);
    
    X_truth(4*i-3:4*i,time_aux)=X_b;
    
    %Sample of target states at the following time steps
    for k=time_aux+1:t_death(i)-1
        delta_tk_i=delta_tk_t(k);
        X_truth(4*i-3:4*i,k)=Euler_simulation(X_truth(4*i-3:4*i,k-1),f_drift,L,q,delta_tk_i,step_size);

    end
    
end

%Evaluation metric parameters
c_gospa=0.1; %Parameter c of the GOSPA metric. We also consider p=2 and alpha=2



%Plot the scenario
figure(5)
clf
hold on
for i=1:N_targets_tot
    plot(X_truth(4*i-3,t_birth(i):t_death(i)-1),X_truth(4*i-1,t_birth(i):t_death(i)-1),'b','Linewidth',1.3)
    plot(X_truth(4*i-3,t_birth(i)),X_truth(4*i-1,t_birth(i)),'ob','Linewidth',1.3,'MarkerFaceColor','b')
    plot(X_truth(4*i-3,t_birth(i):5:t_death(i)-1),X_truth(4*i-1,t_birth(i):5:(t_death(i)-1)),'ob','Linewidth',1.3) 
    
    text(X_truth(4*i-3,t_birth(i))+0,X_truth(4*i-1,t_birth(i))+7,num2str(t_birth(i)))
end
hold off
xlabel('x position (km)')
ylabel('y position (km)')
%axis equal
%axis([0 Area(1) 0 Area(2)])
axis([100 800 6350 6620])
grid on

% figure(6)
% plot(N_targets_alive_t)
% grid on
% xlabel('Time step')
% ylabel('Number of targets')



