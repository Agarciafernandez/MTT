Scenario_number=1;
numtruth = 4; % number of targets

% Initialise model
T = 1;
F = kron(eye(2),[1 T; 0 1]);
Q = 0.01*kron(eye(2),[T^3/3 T^2/2; T^2/2 T]);
H = kron(eye(2),[1 0]);
R =1* eye(2);
chol_R=chol(R)';
p_d=0.90;
p_s=0.99;
Area=[300 300];   
Nsteps=81; %Considered number of time steps in the simulation
l_clutter=10;


intensity_clutter=l_clutter/(Area(1)*Area(2));

%Birth (the birth cover the whole surveillance area)
Ncom_b=1;
weights_b=0.005;
means_b=[100;0;100;0];
P_ini=diag([150 1 150 1].^2);


covs_b(1:4,1:4,1)=P_ini;
%Intensity Poisson prior (time 0)
lambda0=3;


[X_truth,t_birth,t_death]=TrajectoryWilliams15(Scenario_number,Nsteps,F,numtruth,Q,Area);





Nmc=100; %Number of Monte Carlo runs to compute error
c_gospa=10; %Parameter c of the GOSPA metric. We also consider p=2 and alpha=2
gamma_track_metric=1; %Parameter gamma of the metric for sets of trajectories (only when we use this metric)



%Uncomment the following lines to plot this scenario

% figure(5)
% clf
% plot(X_truth(1,t_birth(1):t_death(1)-1),X_truth(3,t_birth(1):t_death(1)-1),'b','Linewidth',1.3)
% hold on
% plot(X_truth(1,t_birth(1)),X_truth(3,t_birth(1)),'xb','Linewidth',1.3)
% plot(X_truth(1,t_birth(1):5:(t_death(1)-1)),X_truth(3,t_birth(1):5:(t_death(1)-1)),'ob','Linewidth',1.3)
% 
% plot(X_truth(5,t_birth(2):t_death(2)-1),X_truth(7,t_birth(2):t_death(2)-1),'r','Linewidth',1.3)
% plot(X_truth(5,t_birth(2)),X_truth(7,t_birth(2)),'xr','Linewidth',1.3)
% plot(X_truth(5,t_birth(2):5:t_death(2)-1),X_truth(7,t_birth(2):5:t_death(2)),'or','Linewidth',1.3)
% 
% 
% plot(X_truth(9,t_birth(3):t_death(3)-1),X_truth(11,t_birth(3):t_death(3)-1),'g','Linewidth',1.3)
% plot(X_truth(9,t_birth(3)),X_truth(11,t_birth(3)),'xg','Linewidth',1.3)
% plot(X_truth(9,t_birth(3):5:t_death(3)-1),X_truth(11,t_birth(3):5:t_death(3)-1),'og','Linewidth',1.3)
% 
% plot(X_truth(13,t_birth(4):t_death(4)-1),X_truth(15,t_birth(4):t_death(4)-1),'black','Linewidth',1.3)
% plot(X_truth(13,t_birth(4)),X_truth(15,t_birth(4)),'xblack','Linewidth',1.3)
% plot(X_truth(13,t_birth(4):5:t_death(4)-1),X_truth(15,t_birth(4):5:t_death(4)-1),'oblack','Linewidth',1.3)
% 
% %axis([0 Area(1) 0 Area(2)])
% hold off
% xlabel('x position (m)')
% ylabel('y position (m)')
% axis equal
% axis([110 180 110 180])
% grid on
% 
