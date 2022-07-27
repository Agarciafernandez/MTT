function [mean_ukf_act,var_ukf_act,z_pred_j,S_j]=linear_kf_update(meank,Pk,A_l,b_l,Omega_l,R,z_real)

z_pred_j=A_l*meank+b_l;
Psi_j=Pk*A_l';

S_j=A_l*Psi_j+R+Omega_l;

Gain=Psi_j/S_j;

mean_ukf_act=meank+Gain*(z_real-z_pred_j);
var_ukf_act=Pk-Gain*Psi_j';

var_ukf_act=(var_ukf_act+var_ukf_act')/2;