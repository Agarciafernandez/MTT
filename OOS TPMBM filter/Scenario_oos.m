rand('seed',9)
randn('seed',9)


%Simulation parameters
Nmc=100; %Number of Monte Carlo runs


%Measurement model parameters

H=kron(eye(2),[1 0]);
R=4*eye(2);
chol_R=chol(R)';
p_d=0.9;
Area=[800 400];
Nsteps=120; %Considered number of time steps in the simulation
l_clutter=10; %Average number of clutter measurements per time step
intensity_clutter=l_clutter/(Area(1)*Area(2)); %Intensity of the clutter

%Parameters of the continous time birth/death process
lambda=0.12;
mu=0.02;

%Mean and covariance matrix at time of appearance
mean_a=[200;3;250;0];
P_a=diag([50 1 50 1].^2);

%We generate the time intervals between different measurements (we take Nsteps measurements)
%delta_tk (which is known but samples from an exponential distribution with parameter mu_delta)
mu_measurement=1; %Note that this is related to time through 1/mu_measurement
delta_tk_t=exprnd(1/mu_measurement,Nsteps,1);
tk=cumsum(delta_tk_t); %Time at which is measurement is obtained


%Process noise intensity (Wiener velocity model)
q=0.2;

%Sampling of the multi-target system model

%We sample the number of new born targets at each time step
lambda_new_born_t=lambda/mu*(1-exp(-mu*delta_tk_t));
N_new_born_t=poissrnd(lambda_new_born_t);
N_targets_tot=sum(N_new_born_t); %Total number of targets in the simulation

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
    
    %Sample of target state at time of birth
    F_k=kron(eye(2),[1 time_lag_i;0 1]);
    Q_k=q*kron(eye(2),[time_lag_i^3/3 time_lag_i^2/2; time_lag_i^2/2 time_lag_i]);
    
    X_b=F_k*X_a+chol(Q_k)'*randn(Nx,1);
    X_truth(4*i-3:4*i,time_aux)=X_b;
    
    %Sample of target states at the following time steps
    for k=time_aux+1:t_death(i)-1
        delta_tk_i=delta_tk_t(k);
        F_k=kron(eye(2),[1 delta_tk_i;0 1]);
        Q_k=q*kron(eye(2),[delta_tk_i^3/3 delta_tk_i^2/2; delta_tk_i^2/2 delta_tk_i]);
        X_truth(4*i-3:4*i,k)=F_k*X_truth(4*i-3:4*i,k-1)+chol(Q_k)'*randn(Nx,1);
    end
    
end

%Evaluation metric parameters
c_gospa=10; %Parameter c of the GOSPA/trajectory metric. We also consider p=2 and alpha=2
gamma_track_metric=1; %Track switching cost for the trajectory metric


%We shuffle the order of the measurements to create out-of-sequence measurements
k_order=1:Nsteps;
N_oos=5;
lambda_oos=1;

for i=1:fix((Nsteps-N_oos)/N_oos)
    change_index=poissrnd(lambda_oos);
    
    %We delay the measurement at this time step
    change_z=k_order(N_oos*i+1:N_oos*i+change_index);
    k_order(N_oos*i+change_index)= k_order(N_oos*i);
    k_order(N_oos*i:N_oos*i+change_index-1)=change_z;
end



%We overwrite tk with the new order
tk=tk(k_order);



%We overwrite delta_tk_t so that it now represents the time difference
%w.r.t current time step
delta_tk_t=zeros(Nsteps,1);
delta_tk_t(1)=tk(1);
max_t=tk(1);

for i=2:Nsteps
    delta_tk_t(i)=tk(i)-max_t;
    if(tk(i)-max_t>0)
        max_t=tk(i);
    end
end




%In-sequence k indices
is_k_indices=delta_tk_t>0;

%Plot time difference between measurements

% figure(1)
% plot(1:Nsteps,delta_tk_t,'b',...
%     1:Nsteps, zeros(Nsteps,1),'--black','Linewidth',1.3)
% grid on
% xlabel('Measurement number')
% ylabel('Time difference')
% axis_to_plot=1:Nsteps;
% hold on
% plot(axis_to_plot(is_k_indices==0),delta_tk_t(is_k_indices==0),'xr','Linewidth',1.3)
% hold off


k_is_series=cumsum(is_k_indices); %Corresponding time step in sequence (to the for loop over k in the main file)
Nsteps_is=sum(is_k_indices); %Number of steps in sequence

tk_is=tk(is_k_indices); %Time corresponding to in-sequence measurements

%In sequence X_truth, births and deaths

X_truth_is=X_truth(:,k_order(is_k_indices));
t_birth_is=zeros(size(t_birth));
t_death_is=zeros(size(t_death));

remove_target_is=[];
for i=1:N_targets_tot
    t_birth_i=find(X_truth_is(4*i-3,:)>0,1);
    
    if(~isempty(t_birth_i))
        
        t_death_i=find(X_truth_is(4*i-3,:)>0,1,'last')+1;
        t_birth_is(i)=t_birth_i;
        t_death_is(i)=t_death_i;
    else
        remove_target_is=[remove_target_is,i];
    end
end

%If one target does not belong to the in-sequence sampling times, it is
%removed.

if(~isempty(remove_target_is))
    t_birth_is(remove_target_is)=[];
    t_death_is(remove_target_is)=[];
    
    X_truth_is(4*remove_target_is-3,:)=[];
    X_truth_is(4*remove_target_is-2,:)=[];
    X_truth_is(4*remove_target_is-1,:)=[];
    X_truth_is(4*remove_target_is,:)=[];
end




%Plot the scenario at the in-sequence times
figure(5)
clf
hold on
for i=1:size(X_truth_is,1)/4
    plot(X_truth_is(4*i-3,t_birth_is(i):t_death_is(i)-1),X_truth_is(4*i-1,t_birth_is(i):t_death_is(i)-1),'b','Linewidth',1.1)
    plot(X_truth_is(4*i-3,t_birth_is(i)),X_truth_is(4*i-1,t_birth_is(i)),'ob','Linewidth',1.1,'MarkerFaceColor','b')
    plot(X_truth_is(4*i-3,t_birth_is(i):5:t_death_is(i)-1),X_truth_is(4*i-1,t_birth_is(i):5:(t_death_is(i)-1)),'ob','Linewidth',1.1)
    
    text(X_truth_is(4*i-3,t_birth_is(i))-20,X_truth_is(4*i-1,t_birth_is(i))+20,num2str(t_birth_is(i)))
end
hold off
xlabel('x position (m)')
ylabel('y position (m)')
axis equal
axis([0 Area(1) 0 Area(2)])
grid on
set(gca,'FontSize',20)


% figure(6)
% plot(N_targets_alive_t)
% grid on
% xlabel('Time step')
% ylabel('Number of targets')
