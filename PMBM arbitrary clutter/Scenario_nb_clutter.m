% Initialise model
T = 1;
F = kron(eye(2),[1 T; 0 1]);
Q = 0.01*kron(eye(2),[T^3/3 T^2/2; T^2/2 T]);
H = kron(eye(2),[1 0]);
R = 4* eye(2);

chol_R=chol(R)';
p_d=0.90;


p_s=0.99;
Area=[300 300];   
Nsteps=81; %Considered number of time steps in the simulation

%Negative binomial clutter parameters
l_clutter=10; %Mean number of clutter measurements
a_nb=20;  %"Dispersion coefficient" a>1


r_nb=l_clutter/(a_nb-1);
p_nb=1/a_nb;



%Birth model
Ncom_b=1;
weights_b=0.1;
means_b=[Area(1)/2;0;Area(2)/2;0];
P_ini=diag([50 1 50 1].^2);


covs_b(1:4,1:4,1)=P_ini;
%Intensity Poisson prior (time 0)
lambda0=5;



%Sampling of the multi-target system model

%We sample the number of new born targets at each time step
lambda_new_born_t=[lambda0;weights_b*ones(Nsteps,1)];
N_new_born_t=poissrnd(lambda_new_born_t);
N_targets_tot=sum(N_new_born_t);


%We go through each target and generate its states, and vectors t_birth and t_death
Nx=4;
Nmc=100; %Number of Monte Carlo runs to compute error


%X_truth matrix with the ground truth set of trajectories/targets
X_truth=zeros(Nx*N_targets_tot,Nsteps);
t_birth=zeros(1,N_targets_tot); %t_birth is the discrete time step at which the target is born
t_death=zeros(1,N_targets_tot); %t_death is the discrete time step the target is already dead

N_new_born_t_prov=N_new_born_t;
time_aux=1; %Auxiliary variable
N_targets_alive_t=zeros(1,Nsteps);

length_trajectories=nbinrnd(1,1-p_s,N_targets_tot,1)+1; %We can sample the length of the targets from a negative binomial distribution. Note we sum one

% length_trajectories=ones(size(length_trajectories));

for i=1:N_targets_tot
    
    while(N_new_born_t_prov(time_aux)==0)
        time_aux=time_aux+1;
    end
    
    N_new_born_t_prov(time_aux)=N_new_born_t_prov(time_aux)-1;
    t_birth(i)=time_aux;
    
    %We need to determine time of death
    
    
    t_death_discrete=min(time_aux+length_trajectories(i),Nsteps+1);
    t_death(i)=t_death_discrete;
    
    N_targets_alive_t(t_birth(i):t_death(i)-1)=N_targets_alive_t(t_birth(i):t_death(i)-1)+1;
    
    %Generation of the samples: Initial time step    
    X_b=means_b+chol(P_ini)'*randn(Nx,1); %This currently considers only one Gaussian component in birth
    X_truth(4*i-3:4*i,time_aux)=X_b;
    
    %Sample of target states at the following time steps
    for k=time_aux+1:t_death(i)-1
        X_prov=F*X_truth(4*i-3:4*i,k-1)+chol(Q)'*randn(Nx,1);
        
        if(X_prov(1)>0 && X_prov(3)>0 && X_prov(1)<Area(1) && X_prov(3)<Area(2)) %We check if the target is in the surveillance area
            X_truth(4*i-3:4*i,k)=X_prov;
        else
            N_targets_alive_t(k:t_death(i)-1)=N_targets_alive_t(k:t_death(i)-1)-1;
            t_death(i)=k;
            break
        end
        
    end
    
end

%Parameter c of the GOSPA metric. We also consider p=2 and alpha=2
c_gospa=10; 



%Plot the scenario

% figure(5)
% clf
% hold on
% for i=1:N_targets_tot
%     plot(X_truth(4*i-3,t_birth(i):t_death(i)-1),X_truth(4*i-1,t_birth(i):t_death(i)-1),'b','Linewidth',1.3)
%     plot(X_truth(4*i-3,t_birth(i)),X_truth(4*i-1,t_birth(i)),'ob','Linewidth',1.3,'MarkerFaceColor','b')
%     plot(X_truth(4*i-3,t_birth(i):10:t_death(i)-1),X_truth(4*i-1,t_birth(i):10:(t_death(i)-1)),'ob','Linewidth',1.3)
%     
%     text(X_truth(4*i-3,t_birth(i))+5,X_truth(4*i-1,t_birth(i)),num2str(t_birth(i)))
%     text(X_truth(4*i-3,t_death(i)-1)+5,X_truth(4*i-1,t_death(i)-1),num2str(t_death(i)-1),'Color','r')
% 
% 
% end
% hold off
% xlabel('x position (m)')
% ylabel('y position (m)')
% axis equal
% axis([0 Area(1) 0 Area(2)])
% grid on



% figure(6)
% plot(N_targets_alive_t)
% grid on
% xlabel('Time step')
% ylabel('Number of targets')
% axis([0 Nsteps 0 14])




% Plot the scenario (x and y components against time) First component x 
% figure(7)
% plotColors = jet(N_targets_tot);
% clf
% hold on
% for i=1:N_targets_tot
%     plot(t_birth(i):t_death(i)-1, X_truth(4*i-3,t_birth(i):t_death(i)-1),'Color',plotColors(i, :),'Linewidth',1.3)
% %    plot(X_truth(4*i-3,t_birth(i):10:t_death(i)-1),X_truth(4*i-1,t_birth(i):10:(t_death(i)-1)),'ob','Linewidth',1.3) 
% end
% axis([0 Nsteps 0 Area(1)])
% xlabel('Time step')
% ylabel('x axis (m)')
% grid on
% pbaspect([2 1 1])

% figure(8)
% clf
% hold on
% for i=1:N_targets_tot
%     plot(t_birth(i):t_death(i)-1, X_truth(4*i-1,t_birth(i):t_death(i)-1),'Color',plotColors(i, :),'Linewidth',1.3) 
%     %plot(X_truth(4*i-3,t_birth(i):10:t_death(i)-1),X_truth(4*i-1,t_birth(i):10:(t_death(i)-1)),'ob','Linewidth',1.3) 
% end
% 
% axis([0 Nsteps 0 Area(2)])
% xlabel('Time step')
% ylabel('y axis (m)')
% grid on
% pbaspect([2 1 1])

% hold off
% xlabel('x position (m)')
% ylabel('y position (m)')
% axis equal
% axis([0 Area(1) 0 Area(2)])
% grid on


% Plot the scenario (range and azimuth components against time) First component x 
% figure(7)
% plotColors = jet(N_targets_tot);
% 
% 
% clf
% hold on
% for i=1:N_targets_tot
%     range_i=sqrt(X_truth(4*i-3,t_birth(i):t_death(i)-1).^2+ X_truth(4*i-1,t_birth(i):t_death(i)-1).^2);
%     plot(t_birth(i):t_death(i)-1, range_i,'Color',plotColors(i, :),'Linewidth',1.3)
% %    plot(X_truth(4*i-3,t_birth(i):10:t_death(i)-1),X_truth(4*i-1,t_birth(i):10:(t_death(i)-1)),'ob','Linewidth',1.3) 
% end
% axis([0 Nsteps 0 400])
% xlabel('Time step')
% ylabel('Range (m)')
% grid on
% pbaspect([2 1 1])
% 
% figure(8)
% clf
% hold on
% for i=1:N_targets_tot
% 
%     azimuth_i=180/pi*atan2(X_truth(4*i-1,t_birth(i):t_death(i)-1),X_truth(4*i-3,t_birth(i):t_death(i)-1));
% 
%     plot(t_birth(i):t_death(i)-1, azimuth_i,'Color',plotColors(i, :),'Linewidth',1.3) 
%     %plot(X_truth(4*i-3,t_birth(i):10:t_death(i)-1),X_truth(4*i-1,t_birth(i):10:(t_death(i)-1)),'ob','Linewidth',1.3) 
% end
% 
% axis([0 Nsteps 0 90])
% xlabel('Time step')
% ylabel('Azimuth (degrees)')
% grid on
% pbaspect([2 1 1])
% hold off
% 
% 
