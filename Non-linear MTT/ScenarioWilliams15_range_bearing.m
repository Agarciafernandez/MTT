Scenario_number=1;
numtruth = 4; % number of targets

% Initialise model
T = 1;
F = kron(eye(2),[1 T; 0 1]);
Q = 0.01*kron(eye(2),[T^3/3 T^2/2; T^2/2 T]);

%Measurement noise models
kappa=1000; %Concentration parameter of VMF distribution


R_range=3;
range_min=10;
range_max=300;
delta_range=range_max-range_min;

p_s=0.99;
Area=[300 300]; %This variable is mainly for display in the range-bearing scenario   
Nsteps=81; %Considered number of time steps in the simulation
l_clutter=10; %Mean number of clutter measurements (Poisson cardinality and uniformly distributed in the surveillance area)

%Calculation of intensity clutter
intensity_clutter=l_clutter/delta_range; %Note besseli(0,0)=1, If the measurement belongs to the field of view the VMF density of clutter has value 1.


%p_d as a function of range
p_d=@(x) exp(-x/(3*range_max));


x_s=[100;100]; %Position of the sensor

%Birth PPP parameters
Ncom_b=1;
weights_b=0.005;


means_b=[180;0;180;0];
P_ini=diag([100 1 100 1].^2);


covs_b(1:4,1:4,1)=P_ini;
%Intensity Poisson prior (time 0)
lambda0=1;



[X_truth,t_birth,t_death]=TrajectoryWilliams15_range_bearing(Scenario_number,Nsteps,F,numtruth,Q,Area);

Nx=4;

c_gospa=10; %Parameter c of the GOSPA metric. We also consider p=2 and alpha=2

Nmc=100;





%Uncomment the following line to show pd_against distance for one target
% figure(6)
% distance_2=sqrt((X_truth(5,t_birth(2):t_death(2)-1)-x_s(1)).^2+(X_truth(7,t_birth(2):t_death(2)-1)-x_s(2)).^2);
% plot(p_d(distance_2))
% xlabel('Time')

%Uncomment the following lines to plot this scenario
% figure(3)
% clf
% [X_p,Y_p]=meshgrid(100:180,100:180);
% distance_p=sqrt((X_p-x_s(1)).^2+(Y_p-x_s(2)).^2);
% 
% pcolor(X_p,Y_p,p_d(distance_p))
% colorbar
% colormap gray
% xlabel('x position (m)')
% ylabel('y position (m)')

% figure (5)
% hold on
% 
% plot(X_truth(1,t_birth(1):t_death(1)-1),X_truth(3,t_birth(1):t_death(1)-1),'b','Linewidth',1.3)
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
% plot(x_s(1),x_s(2),'xr','MarkerSize',10,'Linewidth',3)
% 
% %axis([0 Area(1) 0 Area(2)])
% hold off
% xlabel('x position (m)')
% ylabel('y position (m)')
% axis equal
% axis([99 180 99 180])
% grid on

