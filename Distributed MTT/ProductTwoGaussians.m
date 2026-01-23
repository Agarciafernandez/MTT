function [mean_prod,cov_prod,alpha_prod]=ProductTwoGaussians(mean1,cov1,mean2,cov2)

%This function computes the product of two Gaussians, which is a Gaussian
%multiplied by alpha_prod
%The expression can be found in Eq.(11) in "Williams, J. L. and Maybeck, P.
%S., "Cost-function-based Gaussian mixture reduction for target tracking",
% in Proceedings of the Sixth International Conference of Information Fusion (2003), pp. 1047-1054."

inv_cov1=inv(cov1);
inv_cov2=inv(cov2);

cov_prod=inv(inv_cov1+inv_cov2);

cov_prod=(cov_prod+cov_prod')/2;

mean_prod=cov_prod*(inv_cov1*mean1+inv_cov2*mean2);


alpha_prod=1/sqrt(det(2*pi*(cov1+cov2)))*exp(-0.5*(mean1-mean2)'/(cov1+cov2)*(mean1-mean2));

