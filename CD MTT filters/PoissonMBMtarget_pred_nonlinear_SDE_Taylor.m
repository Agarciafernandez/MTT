function filter_pred=PoissonMBMtarget_pred_nonlinear_SDE_Taylor(filter_upd,f_drift, F_drift,Q_c,p_s,weights_b,means_b,covs_b,delta_t,Nx)

%This function calculates the prediction step with a non-linear SDE

%A. F. García-Fernández, S. Sarkka, "Gaussian multi-target filtering with
%target dynamics driven by a stochastic differential equation", in IEEE
%Transactions on Signal Processing, 2025.


%Author: Ángel F. García-Fernández


%Prediction for the Poisson component
filter_pred.weightPois=p_s*filter_upd.weightPois;
Ncom=length(filter_upd.weightPois);


if(Ncom>0)
    for i=1:Ncom
        mean_act_i=filter_upd.meanPois(:,i);
        cov_act_i=filter_upd.covPois(:,:,i);

        [mean_pred_i,cov_pred_i]=Prediction_SDE_Taylor(mean_act_i,cov_act_i,f_drift, F_drift,Q_c,delta_t,Nx);

        filter_pred.meanPois(:,i)=mean_pred_i;
        filter_pred.covPois(:,:,i)=cov_pred_i;
    end

    %We add PHD of new born targets
    filter_pred.weightPois=cat(2,filter_pred.weightPois,weights_b);
    filter_pred.meanPois=cat(2,filter_pred.meanPois,means_b);
    filter_pred.covPois=cat(3,filter_pred.covPois,covs_b);
else
    filter_pred.weightPois=weights_b;
    filter_pred.meanPois=means_b;
    filter_pred.covPois=covs_b;


end


%Prediction for Bernoulli components
filter_pred.globHyp=filter_upd.globHyp;
filter_pred.globHypWeight=filter_upd.globHypWeight;


Ntracks=length(filter_upd.tracks);

if(Ntracks>0)
    for i=1:Ntracks
        Nhyp_i=length(filter_upd.tracks{i}.eB);
        filter_pred.tracks{i}.t_ini=filter_upd.tracks{i}.t_ini;

        for j=1:Nhyp_i %We go through all hypotheses


            mean_act_ij=filter_upd.tracks{i}.meanB{j};
            cov_act_ij=filter_upd.tracks{i}.covB{j};

            [mean_pred_ij,cov_pred_ij]=Prediction_SDE_Taylor(mean_act_ij,cov_act_ij,f_drift, F_drift,Q_c,delta_t,Nx);

            filter_pred.tracks{i}.meanB{j}=mean_pred_ij;
            filter_pred.tracks{i}.covB{j}=cov_pred_ij;

            filter_pred.tracks{i}.eB(j)=p_s*filter_upd.tracks{i}.eB(j);
            filter_pred.tracks{i}.aHis{j}=filter_upd.tracks{i}.aHis{j};

        end

    end
else

    filter_pred.tracks=cell(0,1);
    filter_pred.globHyp=[];
    filter_pred.globHypWeight=[];

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
