function  [mean_b_k,cov_b_k]=BirthParameters_cd(q,delta_tk,mean_a_p,mean_a_v,P_a_pp,P_a_vv,P_a_pv,mu)

%This function calculate the best Poisson Gaussian fit to the
%continuous-discrete birth model.
%See A F. García-Fernández, S Maskell, "Continuous-discrete multiple target
% ?ltering: PMBM, PHD and CPHD ?lter implementations," in IEEE
% Transactions on Signal Processing, 2020.  DOI: 10.1109/TSP.2020.2968247

%Author: Angel Garcia-Fernandez

%mean_a_s: Mean appearance with the components switched
mean_a_s=[mean_a_p;mean_a_v];
d=length(mean_a_p);
eye_d=eye(d);
zeros_d=zeros(d);

%Moments of the time lag 
E_t=1/mu-delta_tk*exp(-mu*delta_tk)/(1-exp(-mu*delta_tk));
E_t2=1/(1-exp(-mu*delta_tk))*(2/mu^2-exp(-mu*delta_tk)*(delta_tk^2+2*delta_tk/mu+2/mu^2));
E_t3=1/(1-exp(-mu*delta_tk))*(6/mu^3-exp(-mu*delta_tk)*(delta_tk^3+3*delta_tk^2/mu+6*delta_tk/mu^2+6/mu^3));
C_t=E_t2-E_t^2;



%means_b_s: Mean birth with the components switched
mean_b_s=[eye_d,E_t*eye_d;...
    zeros_d,eye_d]*mean_a_s;

P_b_pp=C_t*(mean_a_v*mean_a_v')+q*E_t3/3*eye_d+P_a_pp+E_t*(P_a_pv+P_a_pv')+E_t2*P_a_vv;
P_b_pv=q*E_t2/2*eye_d+P_a_pv+E_t*P_a_vv;
P_b_vv=q*E_t*eye_d+P_a_vv;

%We return the output in the format required by the rest of the code
mean_b_k=zeros(2*d,1);
mean_b_k(1:2:end)=mean_b_s(1:d);
mean_b_k(2:2:end)=mean_b_s(d+1:end);

cov_b_k=zeros(2*d);
cov_b_k(1:2:end,1:2:end)=P_b_pp;
cov_b_k(1:2:end,2:2:end)=P_b_pv;
cov_b_k(2:2:end,1:2:end)=P_b_pv';
cov_b_k(2:2:end,2:2:end)=P_b_vv;
