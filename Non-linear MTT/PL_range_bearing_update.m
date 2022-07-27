function [mean_ukf_act,var_ukf_act,likelihood]=PL_range_bearing_update(meank, Pk,mean_ini,P_ini,z_real,kappa,Ap_kappa,x_s, weights,Nx,Nz,R_kf,Nit_iplf,likelihood_corr_flag,p_d,log_denom_likelihood,threshold_iplf,log_denom_Gauss)

%Update using the IPLF with a VMF bearings measurement and Gaussian range

%If likelihood_corr_flag=1, we calculate the likelihood with the sigma
%points and the correction. Otherwise, we return the standard Gaussian
%approach with constant p_d(evaluated at the prior mean)

%Author: Angel Garcia-Fernandez

W0=weights(1);

flag_convergence=0; %To help debugging

meank_j=mean_ini;
Pk_j=P_ini;


for p=1:Nit_iplf
    
    [A_l,b_l,Omega_l]=SLR_measurement_range_bearings_VM_joint(meank_j,Pk_j,weights,W0,Nx,Nz,kappa,Ap_kappa,x_s);
    
    %KF update (additive noise given by R_kf, only range)
    [mean_ukf_act,var_ukf_act,z_pred,S_pred]=linear_kf_update(meank,Pk,A_l,b_l,Omega_l,R_kf,z_real);
    
    dist_k=dist_kullback(meank_j,Pk_j,mean_ukf_act,var_ukf_act);
    if(dist_k<threshold_iplf)
        flag_convergence=1;
        break;
    end
    
    meank_j=mean_ukf_act;
    Pk_j=var_ukf_act;
    
end


N_points=length(weights);
if(likelihood_corr_flag)
    %We use the sigma point approach w.r.t. the posterior
    %We use log scale to perform these calculations
    %Sigma-point generation
    chol_var_mult=chol((Nx/(1-W0)*var_ukf_act));
    sigma_points=[zeros(Nx,1),chol_var_mult',-chol_var_mult'];
    sigma_points=repmat(mean_ukf_act,1,length(weights))+sigma_points;
    
    %Sigma points p_d
    sigma_points_dist=sqrt((sigma_points(1,:)-x_s(1)).^2+(sigma_points(3,:)-x_s(2)).^2);
    sigma_points_pd_log=log(p_d(sigma_points_dist));
    
    %We first consider VMF part
    sigma_points_pos_norm=zeros(2,N_points);
    sigma_points_pos_norm(1,:)=(sigma_points(1,:)-x_s(1))./sigma_points_dist;
    sigma_points_pos_norm(2,:)=(sigma_points(3,:)-x_s(2))./sigma_points_dist;
    
    %Cross-product with z
    sigma_points_pos_norm_z=sigma_points_pos_norm(1,:)*z_real(1)+sigma_points_pos_norm(2,:)*z_real(2);
    
    %Sigma points true likelihood (log
    sigma_points_l_log=-log_denom_likelihood+kappa*sigma_points_pos_norm_z-(z_real(3)-sigma_points_dist).^2/(2*R_kf(3,3));
    
    %Likelihood_kf_log
    
    Vs= chol(S_pred);
    log_det_S_pred= 2*log(prod(diag(Vs)));
    inv_sqrt_S= inv(Vs);
    inv_S_pred= inv_sqrt_S*inv_sqrt_S';
    quad=-0.5*(z_real-z_pred)'*inv_S_pred*(z_real-z_pred);
    likelihood_kf_log=quad-1/2*log_det_S_pred-log_denom_Gauss;
    
    %Sigma points l_tilde log
    Vs= chol(Omega_l+R_kf);
    log_det_S_pred= 2*log(prod(diag(Vs)));
    inv_sqrt_S= inv(Vs);
    inv_S_pred= inv_sqrt_S*inv_sqrt_S';
    
    sigma_points_l_tilde_log=zeros(1,N_points);
    for i=1:N_points
        z_pred=A_l*sigma_points(:,i)+b_l;
        quad=-0.5*(z_real-z_pred)'*inv_S_pred*(z_real-z_pred);
        sigma_points_l_tilde_log(i)=quad-1/2*log_det_S_pred-log_denom_Gauss;
    end
    
    sigma_points_tot_log=likelihood_kf_log+sigma_points_l_log+sigma_points_pd_log-sigma_points_l_tilde_log;
    
    likelihood=exp(sigma_points_tot_log)*weights';
    
    
else
    %We use the standard approach for calculating the likelihood (with pd evaluated at the prior mean)
    
    %Likelihood_kf_log
    Vs= chol(S_pred);
    log_det_S_pred= 2*log(prod(diag(Vs)));
    inv_sqrt_S= inv(Vs);
    inv_S_pred= inv_sqrt_S*inv_sqrt_S';
    quad=-0.5*(z_real-z_pred)'*inv_S_pred*(z_real-z_pred);
    likelihood_kf_log=quad-1/2*log_det_S_pred-log_denom_Gauss;
    
    %distance at the prior mean
    mean_dist=sqrt((meank(1)-x_s(1)).^2+(meank(3)-x_s(2)).^2);
    
    likelihood=exp(likelihood_kf_log)*p_d(mean_dist);
    
    
end


