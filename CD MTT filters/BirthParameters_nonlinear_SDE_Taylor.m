function  [mean_b,cov_b]=BirthParameters_nonlinear_SDE_Taylor(mean_a,P_a,f_drift, F_drift,Q_c,mu, delta_t)


%This function calculates the best Poisson Gaussian fit to the
%continuous-discrete birth model for a nonlinear SDE using first-order
%Taylor series to approximate the drift function

%A. F. García-Fernández, S. Sarkka, "Gaussian multi-target filtering with
%target dynamics driven by a stochastic differential equation", in IEEE
%Transactions on Signal Processing, 2025.


%Author: Ángel F. García-Fernández

Nx=length(mean_a);


x0=pack(mean_a,P_a,zeros(4,1),zeros(4));


tspan=[0,delta_t];

[t,x_t] = ode45(@(t,x)  odemeancov(t,x,f_drift,F_drift,Q_c,mu,Nx), tspan,x0);

x_final=x_t(end,:)';

[~,~,x_bar_final,Sigma_final]=unpack(x_final,Nx);

%We pre-calculate this exponential factor
exp_factor=mu/(1-exp(-mu*delta_t));

mean_b=exp_factor*x_bar_final;
cov_b=exp_factor*Sigma_final-mean_b*mean_b';

cov_b=(cov_b+cov_b')/2;


end

function x=pack(m,P,x_bar,Sigma)
x=[m;P(:);x_bar;Sigma(:)];
end


function [mean,P,x_bar,Sigma]=unpack(x,Nx)
Nx2=Nx*Nx;

mean=x(1:Nx);
P_list=x(Nx+1:Nx+Nx2);
P=reshape(P_list,Nx,Nx);
x_bar=x(Nx+Nx2+1:Nx2+2*Nx);
Sigma_list=x(Nx2+2*Nx+1:end);

Sigma=reshape(Sigma_list,Nx,Nx);


end

function dxdt = odemeancov(t,x,f_drift,F_drift,Q_c,mu,Nx)

[m,P,x_bar,Sigma]=unpack(x,Nx);

%Differential of the mean
dm=f_drift(m);
%Differential of the covariance
F_grad=F_drift(m);
dP=P*F_grad'+F_grad*P+Q_c;
%Differential of x_bar
dx_bar=m*exp(-mu*t);
%Differential of Sigma
dSigma=(P+m*m')*exp(-mu*t);

dxdt=pack(dm,dP,dx_bar,dSigma);

end

