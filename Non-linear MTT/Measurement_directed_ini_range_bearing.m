function [mean_ini,P_ini]=Measurement_directed_ini_range_bearing(x_s,z_m,kappa,W0)

%We calculated the mean and covariance using sigma points through the
%inverse, see 

% A. F García-Fernández,J. Ralph, P. Horridge, S. Maskell, 
% "A Gaussian filtering method for multi-target tracking with
% nonlinear/non-Gaussian measurements"
% IEEE Transactions on Aerospace and Electronic Systems, 2021, 57, 3539-3548

%Author: Angel Garcia-Fernandez

z_angle=atan2(z_m(2),z_m(1));
z_range=z_m(3);
R=1/kappa;

Nz_angle=1;
Wn=(1-W0)/(2*Nz_angle);
weights_sp=[W0,Wn*ones(1,2*Nz_angle)]; %Weights sigma-points
%Sigma-point generation in angular space
chol_var_mult=chol((Nz_angle/(1-W0)*R));
sigma_points_z=[zeros(Nz_angle,1),chol_var_mult',-chol_var_mult'];
sigma_points_z=repmat(z_angle,1,length(weights_sp))+sigma_points_z;

sigma_points_x=[x_s(1)+z_range*cos(sigma_points_z);x_s(2)+z_range*sin(sigma_points_z)];

%Mean_ini (position elements)
mean_ini_p=sigma_points_x*weights_sp';

P_ini_p=zeros(2);
for i=1:length(weights_sp)
   sub=sigma_points_x(:,i)-mean_ini_p;
   P_ini_p=P_ini_p+weights_sp(i)*(sub*sub'); 
end

%The values of mean_ini(2) and mean_ini(4) do not affect the
%linearisation
mean_ini=zeros(4,1);
mean_ini(1)=mean_ini_p(1);
mean_ini(3)=mean_ini_p(2);

%The values of P_ini(2) and P_ini(4) do not affect the linearisation due to
%the structure
P_ini=zeros(4);
P_ini([1 3],[1 3])=P_ini_p;
P_ini(2,2)=0.1; 
P_ini(4,4)=0.1;