function pd_mean=integral_pd(meank,Pk,p_d,weights,Nx,x_s)

%This function calculates the integral of p_d(x) over a Gaussian (e.g., the
%average p_d)
%Author: Angel F. Garcia-Fernandez
W0=weights(1);
%Sigma-point generation
chol_var_mult=chol((Nx/(1-W0)*Pk));
sigma_points=[zeros(Nx,1),chol_var_mult',-chol_var_mult'];
sigma_points=repmat(meank,1,length(weights))+sigma_points;

sigma_points_dist=sqrt((sigma_points(1,:)-x_s(1)).^2+(sigma_points(3,:)-x_s(2)).^2);
sigma_points_p_d=p_d(sigma_points_dist);

pd_mean=sigma_points_p_d*weights';