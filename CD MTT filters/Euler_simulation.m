function x_k=Euler_simulation(x_ini,f_drift,L,q,time_interval,step_size)


%Euler-Maruyama simulation to obtain a sample from the ODE after a time-interval

Nsteps=ceil(time_interval/step_size); %Note that the last time step has a different time length

x_k_1=x_ini;

cov_noise=q*eye(2)*step_size;
chol_cov_noise=chol(cov_noise)';

x_k_series=zeros(4,Nsteps);


for k=1:Nsteps
    if(k<Nsteps)
        x_k=x_k_1+f_drift(x_k_1)*step_size+L*chol_cov_noise*randn(2,1);
    elseif(k==Nsteps)
        step_size_final=time_interval-(Nsteps-1)*step_size;
        cov_noise=q*eye(2)*step_size_final;
        chol_cov_noise=chol(cov_noise)';
        x_k=x_k_1+f_drift(x_k_1)*step_size_final+L*chol_cov_noise*randn(2,1);
    end
    x_k_1=x_k;

    x_k_series(:,k)=x_k;
end

% figure(1)
% plot(x_k_series(1,1),x_k_series(3,1),'ob')
% hold on
% plot(x_k_series(1,:),x_k_series(3,:))
% hold off
% grid on



