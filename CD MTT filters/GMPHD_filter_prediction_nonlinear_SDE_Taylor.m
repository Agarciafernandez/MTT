function   [weights_k, means_k,covs_k,Ncom_k,logical_active_k]=GMPHD_filter_prediction_nonlinear_SDE_Taylor(weights_u, means_u,covs_u,Ncom_u,logical_active_u,f_drift, F_drift,Q_c,p_s,delta_t,Nx)

%Prediction for the CD-GMPHD filter with a non-linear SDE as the dynamic
%model

%Author: Ángel F. García-Fernández

Ncom_k=Ncom_u;
logical_active_k=logical_active_u;

index_active_u=find(logical_active_u==1);

weights_k=zeros(size(weights_u));
means_k=zeros(size(means_u));
covs_k=zeros(size(covs_u));

for i=1:Ncom_u
    mean_i=means_u(:,index_active_u(i));
    cov_i=covs_u(:,:,index_active_u(i));
    weight_i=weights_u(index_active_u(i));

    [mean_i_pred,cov_i_pred]=Prediction_SDE_Taylor(mean_i,cov_i,f_drift, F_drift,Q_c,delta_t,Nx);

    weight_i_pred=weight_i*p_s;

    weights_k(index_active_u(i))=weight_i_pred;
    means_k(:,index_active_u(i))=mean_i_pred;
    covs_k(:,:,index_active_u(i))=cov_i_pred;

end

end

function [mean_pred,cov_pred]=Prediction_SDE_Taylor(mean_act,cov_act,f_drift, F_drift,Q_c,delta_t,Nx)

%This function provides the prediction step for a Gaussian density
x0=pack(mean_act,cov_act);

tspan=[0,delta_t];

[t,x_t] = ode45(@(t,x)  odemeancov(t,x,f_drift,F_drift,Q_c,Nx), tspan,x0);

x_final=x_t(end,:)';

[mean_pred,cov_pred]=unpack(x_final,Nx);


end

function x=pack(m,P)
x=[m;P(:)];
end


function [mean,P]=unpack(x,Nx)
Nx2=Nx*Nx;

mean=x(1:Nx);
P_list=x(Nx+1:Nx+Nx2);
P=reshape(P_list,Nx,Nx);
end


function dxdt = odemeancov(t,x,f_drift,F_drift,Q_c,Nx)

[m,P]=unpack(x,Nx);

%Differential of the mean
dm=f_drift(m);
%Differential of the covariance
F_grad=F_drift(m);
dP=P*F_grad'+F_grad*P+Q_c;


dxdt=pack(dm,dP);

end

