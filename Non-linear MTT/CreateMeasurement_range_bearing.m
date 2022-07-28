function z=CreateMeasurement_range_bearing(X_truth_k,t_birth,t_desth,p_d,l_clutter,k,Nx,x_s,kappa,R_range,range_min,delta_range)
%This function generates range bearings measurements where range is
%distributed as Gaussian and bearings as VMF.

%Author: Angel F. Garcia-Fernandez
X_truth_k_r=reshape(X_truth_k,Nx,size(X_truth_k,2)*size(X_truth_k,1)/Nx);

range_targets=sqrt((X_truth_k_r(1,:)-x_s(1)).^2+(X_truth_k_r(3,:)-x_s(2)).^2);

p_d_targets=p_d(range_targets);

%Targets must be alive to create a measurement
index_targets=and(t_birth<=k,t_desth>k);

detected_targets=rand(1,length(index_targets))<p_d_targets;
detected_targets=and(detected_targets,index_targets);

N_detected_targets=sum(detected_targets);
N_clutter_measurements=poissrnd(l_clutter);

%We create the range bearing measurements (We first put bearings,
%indicated as a 2-D measurement and then range)

%Angular measurement
mean_z=[X_truth_k_r(1,detected_targets)-x_s(1);X_truth_k_r(3,detected_targets)-x_s(2)];
mean_theta=atan2(mean_z(2,:),mean_z(1,:));
%We generate VMF
z_real_angle=zeros(1,N_detected_targets);
for i=1:N_detected_targets
    z_real_angle(i) = circ_vmrnd(mean_theta(i), kappa, 1);
end
%We generate range
mean_z_range=range_targets(detected_targets);
z_real_range=mean_z_range+sqrt(R_range)*randn(1,N_detected_targets);

%We generate teh 3-D measurement (with direction and range)
z_targets=[cos(z_real_angle);sin(z_real_angle);z_real_range];

%We generate clutter
angle_uniform=2*pi*rand(1,N_clutter_measurements);
range_uniform=range_min+delta_range*rand(1,N_clutter_measurements);

z_clutter=[cos(angle_uniform);sin(angle_uniform);range_uniform];

z=[z_targets,z_clutter];





