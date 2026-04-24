function beta_score=CalibrateBetaScore(c_gospa,p,rho_qmetric,n_fa,score_n_fa)

%Function to calibrate the beta parameter for the score function based on
%the GOSPA and the T-GOSPA q-metrics

%n_fa: number of false targets
%score_n_fa: score for this number of false targets

%This is implemented for the sigmoid-based metric-preserving mapping

%A. F. García-Fernández, J. Gu, L. Svensson, Y. Xia, J. Krejcí, O. Kost, O. Straka “GOSPA and T-GOSPA quasi-metrics for evaluation of multi-object tracking algorithms,” 
% accepted in IEEE Transactions on Aerospace and Electronic Systems, 2026.


Normalised_metric_value=inverse_sigmoid_mapping(1-score_n_fa);


Metric_value=c_gospa*nthroot(rho_qmetric*n_fa,p);

beta_score=Metric_value/Normalised_metric_value;




end




% 
% function y=sigmoid_mapping(x)
% %Sigmoid metric preserving mapping
% 
% y=2/(1+exp(-x))-1;
% 
% end



function x=inverse_sigmoid_mapping(y)
%Sigmoid metric preserving mapping

x=-log(2/(y+1)-1);

end