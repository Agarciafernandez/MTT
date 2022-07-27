function [z_pred_ukf,var_pred_ukf]=Predicted_range_bearings_VM(meank,Pk,weights,W0,Nx,Nz,kappa,Ap_kappa,x_s)

%This function calculates the expected mean and covariance of a von Mises
%Fisher distribution (also with range), without accounting for the range
%noise.

%Author: Angel Garcia-Fernandez


Nz_angle=Nz-1;

%Sigma-point generation
chol_var_mult=chol((Nx/(1-W0)*Pk));
sigma_points=[zeros(Nx,1),chol_var_mult',-chol_var_mult'];
sigma_points=repmat(meank,1,length(weights))+sigma_points;


%First the angle (transformation h) and then range (transformation r)
%Transformation through h
sigma_points_r=sqrt((sigma_points(1,:)-x_s(1)).^2+(sigma_points(3,:)-x_s(2)).^2);
sigma_points_h=[sigma_points(1,:)-x_s(1);sigma_points(3,:)-x_s(2)]./repmat(sigma_points_r,2,1);
mean_h=sigma_points_h*weights';

%Range
mean_r=sigma_points_r*weights';
sub_z=sigma_points_r-mean_r;
cov_r=sum(weights.*sub_z.^2);



mean_hh=zeros(Nz_angle);%E[h*h']
cov_hr=zeros(Nz_angle,1);
for j=1:length(weights)
    mean_hh=mean_hh+weights(j)*(sigma_points_h(:,j)*sigma_points_h(:,j)');
    cov_hr=cov_hr+weights(j)*(sigma_points_h(:,j)-mean_h)*(sigma_points_r(j)-mean_r);
end

%Now we can compute the moments of the transformed variable (using the VM
%adjustments to transfrom from h to measurement)
p=2; %p=2 as we are on 2D
z_pred_ukf=[Ap_kappa*mean_h;mean_r];

cov_angle=(Ap_kappa)^2*(mean_hh-mean_h*mean_h')+...
    +Ap_kappa/kappa*eye(Nz_angle)+(1-(Ap_kappa)^2-p*Ap_kappa/kappa)*mean_hh';
cov_angle_range=Ap_kappa*cov_hr;

var_pred_ukf=[cov_angle,cov_angle_range;cov_angle_range',cov_r];