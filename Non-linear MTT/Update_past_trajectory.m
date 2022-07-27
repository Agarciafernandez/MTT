function [mean_u,P_u]=Update_past_trajectory(mean_i,cov_i,x_u_c,P_u_c,Nx)

%Updated moments of the whole trajectory given the updated moments of the target at the current time step, see (68)-(70) in "Gaussian MAP Filtering
%Using Kalman Optimization" and 

%F. Beutler, M. F. Huber, and U. D. Hanebeck, “Gaussian filtering using
%state decomposition methods," in 12th International Conference on
%Information Fusion, July 2009, pp. 579–586.


%mean_i and cov_i prior mean and covariance
%x_u_c and P_u_c are updated mean and covariance for the current target
%state

%Mean current and past (within window of update)
mean_i_c=mean_i(end-Nx+1:end);
mean_i_p=mean_i(1:end-Nx);

%Covariance current state and past states, and cross-covariance

cov_i_c=cov_i(end-Nx+1:end,end-Nx+1:end);
cov_i_p=cov_i(1:end-Nx,1:end-Nx);
cov_i_cp=cov_i(end-Nx+1:end,1:end-Nx);

%Updated moments, see (68)-(70) in "Gaussian MAP Filtering
%Using Kalman Optimization"
%Gain
G=cov_i_cp'/cov_i_c;
x_u_p=mean_i_p+G*(x_u_c-mean_i_c);
P_u_p=cov_i_p-G*(cov_i_c-P_u_c)*G';
P_u_cp=P_u_c*G';

mean_u=[x_u_p;x_u_c];
P_u=[P_u_p P_u_cp';P_u_cp P_u_c];
P_u=(P_u+P_u')/2;