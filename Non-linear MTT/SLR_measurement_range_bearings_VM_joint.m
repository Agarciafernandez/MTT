function [A_l,b_l,Omega_l]=SLR_measurement_range_bearings_VM_joint(meank,Pk,weights,W0,Nx,Nz,kappa,Ap_kappa,x_s)
%This function computes the statistical linear regression (SLR) of conditional moments for Von Mises distribution.
% It compute the joint linearisation with range and bearings. It is also possible to
%do it indepedently for each type of measurement.

%First we compute the moments of h(x), which is the function that project
%the state to the unit circle

%Author: Angel F. Garcia-Fernandez

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
cov_xr=sum((sigma_points-meank).*repmat(weights.*sub_z,4,1),2);

%Sigma_points with no mean
sigma_points_h_centered=sigma_points_h-mean_h;
sigma_points_r_centered=sigma_points_r-mean_r;

mean_hh=zeros(Nz_angle);%E[h*h']
cov_xh=zeros(Nx,Nz_angle);
cov_hr=zeros(Nz_angle,1);

for j=1:length(weights)
    mean_hh=mean_hh+weights(j)*(sigma_points_h(:,j)*sigma_points_h(:,j)');
    cov_xh=cov_xh+weights(j)*((sigma_points(:,j)-meank)*sigma_points_h_centered(:,j)');
    cov_hr=cov_hr+weights(j)*sigma_points_h_centered(:,j)*sigma_points_r_centered(j);
end

%Now we can compute the moments of the transformed variable (using the VM
%adjustments to transfrom from h to measurement)
p=2; %p=2 as we are on 2D
z_pred_ukf=[Ap_kappa*mean_h;mean_r];

var_xz_ukf=[Ap_kappa*cov_xh,cov_xr];
cov_angle=(Ap_kappa)^2*(mean_hh-mean_h*mean_h')+...
    +Ap_kappa/kappa*eye(Nz_angle)+(1-(Ap_kappa)^2-p*Ap_kappa/kappa)*mean_hh';
cov_angle_range=Ap_kappa*cov_hr;

var_pred_ukf=[cov_angle,cov_angle_range;cov_angle_range',cov_r];

%Statistical linearisaltion
A_l=var_xz_ukf'/Pk;
b_l=z_pred_ukf-A_l*meank;
%Omega_l=var_pred_ukf-A_l*Pk*A_l';

Omega_l=var_pred_ukf-var_xz_ukf'*A_l'; %Faster






