function  [mean_b,cov_b]=BirthParameters_linear_SDE(mean_a,P_a,A,u,Q_c,mu, delta_t)


%This function calculates the best Poisson Gaussian fit to the
%continuous-discrete birth model for a linear SDE in 

%A. F. García-Fernández, S. Sarkka, "Gaussian multi-target filtering with
%target dynamics driven by a stochastic differential equation", in IEEE
%Transactions on Signal Processing, 2025.

%This version computes several terms with a single matrix exponential for a
%faster implementation

%Author: Ángel F. García-Fernández


Nx=length(mean_a);


%We pre-calculate this exponential factor
exp_fact=1/(1-exp(-mu*delta_t));
exp_fact1=mu*exp_fact;

%We calculate the mean at the time of birth

A_b_u=[A-mu*eye(Nx),u,zeros(Nx,1);...
    zeros(1,Nx), -mu , 1;...
    zeros(1,Nx), 0, 0];

exp_A_b_u=expm(A_b_u*delta_t);

mean_b_u1=exp_fact1*[eye(Nx), zeros(Nx,2)]*exp_A_b_u*[zeros(Nx+1,1);1];

mean_b_1=exp_fact1*Calculate_Integral_type1(A-mu*eye(Nx),mean_a,delta_t);

mean_b=mean_b_1+mean_b_u1;

%We calculate the covariance at the time of birth

%We first calculate E_C_x_1 (first term)
E_C_x_1=-exp(-mu*delta_t)*exp_fact*Calculate_Integral_type2(A,Q_c,delta_t);

%Now we combina E_C_x_2 (second term with Sigma_xx)
C_b=Q_c+mu*P_a;

E_C_x_2_Sigma_xx=Calculate_Integral_type2(A-mu/2*eye(Nx),C_b*exp_fact+exp_fact1*(mean_a*mean_a'),delta_t);

%We calculate Sigma_xu+Sigma_uu2 

Sigma_part1a=Calculate_Integral_type3(A-mu*eye(Nx),exp_fact1*mean_a*u'+exp_fact*u*u',A',delta_t);
%We calculate Sigma_xu+ Sigma_uu2+Sigma_xu'+Sigma_uu2' 
Sigma_part1=Sigma_part1a+Sigma_part1a';


Integral_uu1=Calculate_Integral_type1(A,u,delta_t);
Sigma_uu1=Integral_uu1*Integral_uu1'*exp(-mu*delta_t)*exp_fact;




cov_b=E_C_x_1+E_C_x_2_Sigma_xx+Sigma_part1-Sigma_uu1-mean_b*mean_b';

cov_b=(cov_b+cov_b')/2;

