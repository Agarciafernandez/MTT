rand('seed',9)
randn('seed',9)

%Simulation parameters
Nmc=100; %Number of Monte Carlo runs

%Measurement parameters

H=kron(eye(2),[1 0]);
R=4*eye(2);
chol_R=chol(R)';
p_d=0.9;
Area=[600 400];
Nsteps=100; %Considered number of time steps in the simulation
l_clutter=10; %Average number of clutter measurements per time step
intensity_clutter=l_clutter/(Area(1)*Area(2)); %Intensity of the clutter

%Parameters of the dynamic model
p_s=0.99;
T=1;
F=kron(eye(2),[1 T; 0 1]);
Q=0.01*kron(eye(2),[T^3/3 T^2/2; T^2/2 T]);

chol_Q=chol(Q)';

%Birth
Ncom_b=1;
weights_b=0.08;
means_b=[300;3;170;1];
covs_b=diag([160 1 100 1].^2);

%Spawning modes (two spawning modes with same F and Q but an offset parallel to orthogonal to
%direction of movement)
d_spawning=5; %Spawning distance
p_s1=0.01; %Spawning probabilities
p_s2=0.01;


F1=[1 0 0 -T;
    0 0 0 -1;
    0 T 1 0;
    0 1 0 0]; % Change of velocities in x and y and change of sign in x

F2=[1 0 0 T;
    0 0 0 1;
    0 -T 1 0;
    0 -1 0 0];


%Sampling of the multi-target system model
%We go through each target and generate its states, and vectors t_birth and t_death
Nx=4;

X_truth=[];
t_birth=[];
t_death=[];
father_spawning=[]; %If one trajectory has spawned, it stores the index of the target that has spawned

index_new=1; %Index to allocate new targets (Number of targets +1)

for k=1:Nsteps
    if(k>1)
        %We go through previous targets
        for i=1:index_new-1
            
            X_k_1=X_truth(4*i-3:4*i,k-1);
            
            if(sum(X_k_1==0)<Nx)
                %Target is alive at k-1 (otherwise, we do not do anything
                
                %Survival
                aux=rand(3,1); %Auxiliary variable for survival and two spawning modes
                if(aux(1)<p_s)
                    %Target survives
                    X_k=F*X_k_1+chol_Q*randn(Nx,1);
                    X_truth(4*i-3:4*i,k)=X_k;
                else
                    %Target dies, we set t_death
                    t_death(i)=k;
                end
                
                %Unit vector proportional to velocity (for spawning)
                    vel=sqrt(X_k_1(2)^2+X_k_1(4)^2);
                    unit_vector=[X_k_1(2)/vel,0,X_k_1(4)/vel,0]';
                    
                %Spawning mode 1
                if(aux(2)<p_s1)
                    %Target spawns with mode 1                    
                    d_offset=d_spawning*[-unit_vector(3),0,unit_vector(1),0]';
                    X_k=d_offset+F1*X_k_1+chol_Q*randn(Nx,1);
                    
                    X_truth(4*index_new-3:4*index_new,k)=X_k;
                    
                    t_birth(index_new)=k;
                    father_spawning(index_new)=i;
                    index_new=index_new+1;
                end
                
                %Spawning mode 2
                if(aux(3)<p_s2)
                    %Target spawns with mode 2
                    d_offset=d_spawning*[unit_vector(3),0,-unit_vector(1),0]'; %Signed changed w.r.t. spawning mode 1
                    X_k=d_offset+F2*X_k_1+chol_Q*randn(Nx,1);
                    
                    X_truth(4*index_new-3:4*index_new,k)=X_k;
                    
                    t_birth(index_new)=k;
                    father_spawning(index_new)=i;
                    index_new=index_new+1;
                end
                
                
            end
        end
        
    end
    
    
    %New targets
    N_new_targets=poissrnd(weights_b);
    
    for i=1:N_new_targets
        X_truth(4*index_new-3:4*index_new,k)=means_b+chol(covs_b)'*randn(Nx,1);
        t_birth(index_new)=k;
        father_spawning(index_new)=0; %This is the root node
        index_new=index_new+1;
    end
    
end

%We fill out the missing values in t_death
t_death(t_death==0)=Nsteps+1;
t_death=[t_death,(Nsteps+1)*ones(1,length(t_birth)-length(t_death))];



%Evaluation metric parameters
c_gospa=10; %Parameter c of the GOSPA metric. We also consider p=2 and alpha=2
gamma_track_metric=1; %Parameter gamma of the metric for sets of trajectories (only when we use this metric)



N_targets_tot=length(t_birth);

%Plot the scenario
% figure(5)
% clf
% hold on
% for i=1:N_targets_tot
%     plot(X_truth(4*i-3,t_birth(i):t_death(i)-1),X_truth(4*i-1,t_birth(i):t_death(i)-1),'b','Linewidth',1.2)
%     
%     plot(X_truth(4*i-3,t_birth(i):10:t_death(i)-1),X_truth(4*i-1,t_birth(i):10:(t_death(i)-1)),'ob','Linewidth',1.2,'MarkerSize',3)
%     
%     if(father_spawning(i)==0)
%         plot(X_truth(4*i-3,t_birth(i)),X_truth(4*i-1,t_birth(i)),'ob','Linewidth',1.2,'MarkerFaceColor','b','MarkerSize',3)
%     else
%         plot(X_truth(4*i-3,t_birth(i)),X_truth(4*i-1,t_birth(i)),'sr','Linewidth',1.2,'MarkerFaceColor','r','MarkerSize',4)     
%     end
%     
%     text(X_truth(4*i-3,t_birth(i))-10,X_truth(4*i-1,t_birth(i))+10,num2str(t_birth(i)),'FontSize',8)
% end
% 
% %Plot links for spawning parent-child
% for i=1:N_targets_tot   
%     if(father_spawning(i)>0)
%         link_x=[X_truth(4*father_spawning(i)-3,t_birth(i)-1);...
%             X_truth(4*i-3,t_birth(i))];
%          link_y=[X_truth(4*father_spawning(i)-1,t_birth(i)-1);...
%             X_truth(4*i-1,t_birth(i))];
%         plot(link_x,link_y,'b','Linewidth',1.2)
%     end
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